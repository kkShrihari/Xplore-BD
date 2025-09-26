"""
gene.py — Unified gene information bridge across NCBI, Ensembl, GTEx, ENCODE, DisGeNET, etc.
"""

import os
import requests
import pandas as pd
from dataclasses import dataclass
import xml.etree.ElementTree as ET
from typing import List, Any, Dict, Optional


@dataclass
class Config:
    timeout: float = 30.0
    request_delay: float = 0.2
    ensembl_base: str = "https://rest.ensembl.org"
    ncbi_base: str = "https://api.ncbi.nlm.nih.gov"
    gtex_base: str = "https://gtexportal.org/api/v2"
    encode_base: str = "https://www.encodeproject.org"
    disgenet_base: str = "https://www.disgenet.org/api"
    disgenet_token: Optional[str] = os.getenv("DISGENET_API_KEY")  # read from env


class GeneBridge:
    """Bridge for gene annotation, expression, variants, regulation, disease associations."""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "Explore_BioDB/0.1",
            "Accept": "application/json"
        })

    # ---------- Helpers ----------
    def _resolve_symbol(self, symbol: str, species: str = "human") -> Optional[str]:
        """Resolve symbol to versioned Ensembl ID (e.g. ENSG00000012048.15)."""
        url = "https://mygene.info/v3/query"
        params = {
            "q": f"symbol:{symbol}",
            "species": species,
            "fields": "ensembl.gene,ensembl.transcript"
        }
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        if not r.ok:
            return None
        hits = r.json().get("hits", [])
        if not hits:
            return None
        ensembl_info = hits[0].get("ensembl")
        if isinstance(ensembl_info, list):
            return ensembl_info[0].get("gene")
        elif isinstance(ensembl_info, dict):
            return ensembl_info.get("gene")
        return None

    
    def _resolve_entrez_id(self, gene: str) -> str:
        """Resolve gene symbol to Entrez Gene ID using NCBI E-utilities."""
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {"db": "gene", "term": f"{gene}[sym] AND human[orgn]", "retmode": "json"}
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        if not ids:
            raise ValueError(f"Could not resolve {gene} to Entrez Gene ID")
        return ids[0]

    # ---------- Gene Annotation ----------
    def get_annotation(self, gene: str, source: str = "ensembl") -> Dict[str, Any]:
        """Fetch gene annotation, function, GO terms, orthologs."""
        if source == "ensembl":
            ens_id = self._resolve_symbol(gene) if not gene.startswith("ENSG") else gene
            if not ens_id:
                raise ValueError(f"Could not resolve symbol {gene} to Ensembl ID")
            url = f"{self.config.ensembl_base}/lookup/id/{ens_id}"
            params = {"expand": "0"}
            r = self.session.get(url, params=params, timeout=self.config.timeout)
            r.raise_for_status()
            return r.json()
        elif source == "ncbi":
            entrez_id = self._resolve_entrez_id(gene)
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            params = {"db": "gene", "id": entrez_id, "retmode": "json"}
            r = self.session.get(url, params=params, timeout=self.config.timeout)
            r.raise_for_status()
            data = r.json()
            record = data["result"][entrez_id]

            # 🔎 Keep only the most recent GRCh38 location
            for loc in record.get("locationhist", []):
                if loc.get("assemblyaccver", "").startswith("GCF_000001405"):  # GRCh38
                    record["genomic_location"] = {
                        "assembly": loc.get("assemblyaccver"),
                        "chr_accession": loc.get("chraccver"),
                        "start": loc.get("chrstart"),
                        "end": loc.get("chrstop"),
                    }
                    break

            record.pop("locationhist", None)
            return record
        else:
            raise ValueError("Unsupported source")

    
    # ---------- Gene Expression ----------
    def get_expression(self, gene: str, organism: str = "Homo sapiens") -> Dict[str, Any]:
        import pandas as pd

        ens_id = self._resolve_symbol(gene, organism)
        if not ens_id:
            raise ValueError(f"Could not resolve {gene} to Ensembl ID")

        result = {
            "ensembl_id": ens_id,
            "symbol": gene,
            "organism": organism,
            "expression": [],
            "count_matrix": "https://gtexportal.org/home/datasets",
            "fastq_links": []
        }

        # ---------- Step 1: Try GTEx API ----------
        url = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
        try:
            params = {"gencodeId": ens_id, "datasetId": "gtex_v8"}
            r = self.session.get(url, params=params, timeout=self.config.timeout)
            r.raise_for_status()
            data = r.json()
            expr_data = data.get("medianGeneExpression", [])
        except Exception:
            expr_data = []

        # ---------- Step 2: Fallback to GCT file ----------
        if not expr_data:
            try:
                gct_url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
                df = pd.read_csv(gct_url, sep="\t", skiprows=2)
                row = df[(df["Description"] == gene) | (df["Name"].str.split(".").str[0] == ens_id)]
                if not row.empty:
                    row = row.iloc[0]
                    expr_data = [
                        {"tissue": col, "value": float(row[col]), "unit": "TPM"}
                        for col in df.columns[2:]
                    ]
            except Exception as e:
                print(f"Warning: Could not fetch GTEx matrix fallback ({e})")

        # ---------- Step 3: Save expression ----------
        if expr_data:
            result["expression"] = expr_data
        else:
            print(f"Warning: No expression data found for {gene}")

        # ---------- Step 4: Add FASTQ links (ENA) ----------
        try:
            ena_url = "https://www.ebi.ac.uk/ena/portal/api/search"
            ena_params = {
                "result": "read_run",
                "query": "study_accession=PRJNA758389",
                "fields": "run_accession,fastq_ftp",
                "format": "json"
            }
            ena_r = self.session.get(ena_url, params=ena_params, timeout=self.config.timeout)
            if ena_r.ok:
                runs = ena_r.json()
                result["fastq_links"] = [
                    f["fastq_ftp"] for f in runs if "fastq_ftp" in f
                ][:5]
        except Exception as e:
            print(f"Warning: Could not fetch FASTQ links ({e})")

        return result



    # ---------- Gene Variants ----------
    def get_variants(
        self,
        gene: str,
        source: str = "clinvar",
        min_results: int = 3
    ) -> List[Dict[str, Any]]:
        """
        Fetch gene variants from ClinVar (default) or dbSNP.
        - min_results: number of variants to return (default=3, max=20).
        Uses efetch (VCV XML) for detailed data, with esummary fallback (ClinVar).
        """

        entrez_id = self._resolve_entrez_id(gene)
        results = []
        limit = min(20, max(1, min_results))

        def normalize(value):
            """Convert None/Unknown/empty to 'not specified'."""
            if value in (None, "", "Unknown"):
                return "not specified"
            return value

        # ---------- ClinVar ----------
        if source == "clinvar":
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {"db": "clinvar", "term": f"{entrez_id}[geneid]", "retmode": "json", "retmax": 50}
            r = self.session.get(esearch_url, params=params, timeout=self.config.timeout)
            r.raise_for_status()
            idlist = r.json().get("esearchresult", {}).get("idlist", [])

            if not idlist:
                return []

            for vid in idlist:
                if len(results) >= limit:
                    break

                # Try efetch VCV XML
                efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                params = {"db": "clinvar", "id": vid, "rettype": "vcv", "retmode": "xml"}
                r2 = self.session.get(efetch_url, params=params, timeout=self.config.timeout)

                if r2.ok and "<VariationArchive" in r2.text:
                    try:
                        root = ET.fromstring(r2.text)
                        va = root.find(".//VariationArchive")
                        if va is None:
                            continue

                        record = {
                            "variation_id": normalize(va.attrib.get("VariationID")),
                            "accession": normalize(va.attrib.get("Accession")),
                            "title": normalize(va.attrib.get("VariationName")),
                            "variant_type": normalize(va.attrib.get("VariationType")),
                            "cytogenetic_location": "not specified",
                            "genomic_locations": {},
                            "hgvs": [],
                            "protein_change": "not specified",
                            "molecular_consequence": [],
                            "clinical_significance": "not specified",
                            "review_status": "not specified",
                            "last_evaluated": "not specified",
                            "conditions": [],
                            "submissions": [],
                            "date_created": normalize(va.attrib.get("DateCreated")),
                            "last_updated": normalize(va.attrib.get("DateLastUpdated")),
                        }

                        # Cytogenetic location
                        cytoloc = va.find(".//CytogeneticLocation")
                        if cytoloc is not None and cytoloc.text:
                            record["cytogenetic_location"] = cytoloc.text

                        # Sequence locations
                        for loc in va.findall(".//SequenceLocation"):
                            asm = loc.attrib.get("Assembly")
                            acc = loc.attrib.get("Accession")
                            start = loc.attrib.get("start")
                            stop = loc.attrib.get("stop")
                            if asm and acc and start and stop:
                                record["genomic_locations"][asm] = f"{acc}:{start}-{stop}"

                        # HGVS & protein
                        for hgvs in va.findall(".//HGVS//Expression"):
                            if hgvs.text:
                                record["hgvs"].append(hgvs.text)
                        for prot in va.findall(".//ProteinExpression/Expression"):
                            if prot.text:
                                record["protein_change"] = prot.text

                        # Molecular consequences
                        for mc in va.findall(".//MolecularConsequence"):
                            if "Type" in mc.attrib:
                                record["molecular_consequence"].append(mc.attrib["Type"])

                        # Clinical info
                        cs = va.find(".//Classifications//GermlineClassification//Description")
                        if cs is not None and cs.text:
                            record["clinical_significance"] = cs.text
                        rs = va.find(".//Classifications//GermlineClassification//ReviewStatus")
                        if rs is not None and rs.text:
                            record["review_status"] = rs.text
                        le = va.find(".//Classifications//GermlineClassification")
                        if le is not None:
                            record["last_evaluated"] = normalize(le.attrib.get("DateLastEvaluated"))

                        # Conditions
                        for cond in va.findall(".//ClassifiedConditionList//ClassifiedCondition"):
                            if cond.text:
                                record["conditions"].append(cond.text)

                        # Submissions
                        for scv in va.findall(".//ClinVarAccession[@Type='SCV']"):
                            record["submissions"].append(scv.attrib.get("Accession"))

                        results.append(record)
                        continue
                    except Exception:
                        pass

                # fallback esummary
                summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
                params = {"db": "clinvar", "id": vid, "retmode": "json"}
                r3 = self.session.get(summary_url, params=params, timeout=self.config.timeout)
                if not r3.ok:
                    continue

                record = r3.json().get("result", {}).get(str(vid), {})
                if not record:
                    continue

                germline = record.get("germline_classification", {})

                results.append({
                    "variation_id": normalize(vid),
                    "accession": normalize(record.get("accession")),
                    "title": normalize(record.get("title")),
                    "variant_type": normalize(record.get("obj_type")),
                    "cytogenetic_location": normalize(record.get("band")),
                    "genomic_locations": {"chr": normalize(record.get("chr_sort"))},
                    "hgvs": [],
                    "protein_change": normalize(record.get("protein_change")),
                    "molecular_consequence": record.get("molecular_consequence_list", []),
                    "clinical_significance": normalize(germline.get("description")),
                    "review_status": normalize(germline.get("review_status")),
                    "last_evaluated": normalize(germline.get("last_evaluated")),
                    "conditions": [
                        t.get("trait_name") for t in germline.get("trait_set", []) if "trait_name" in t
                    ] or [],
                    "submissions": record.get("supporting_submissions", {}).get("scv", []) or [],
                    "date_created": normalize(record.get("first_in_clinvar")),
                    "last_updated": normalize(record.get("last_updated")),
                })

        # ---------- dbSNP ----------
        elif source == "dbsnp":
            esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                "db": "snp",
                "term": f"{entrez_id}[GeneID]",
                "retmode": "json",
                "retmax": 50
            }
            r = self.session.get(esearch_url, params=params, timeout=self.config.timeout)
            r.raise_for_status()
            idlist = r.json().get("esearchresult", {}).get("idlist", [])

            if not idlist:
                return []

            for rsid in idlist[:limit]:
                efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                params = {"db": "snp", "id": rsid, "rettype": "docsum", "retmode": "xml"}
                r2 = self.session.get(efetch_url, params=params, timeout=self.config.timeout)
                if not r2.ok:
                    continue

                try:
                    root = ET.fromstring(r2.text)
                    docsum = root.find(".//DocumentSummary")
                    if docsum is None:
                        continue

                    record = {
                        "variation_id": normalize(docsum.attrib.get("uid")),
                        "accession": f"rs{normalize(docsum.findtext('SNP_ID'))}",
                        "title": normalize(docsum.findtext("DOCSUM")),
                        "variant_type": normalize(docsum.findtext("SNP_CLASS")),
                        "alleles": [],
                        "clinical_significance": normalize(docsum.findtext("CLINICAL_SIGNIFICANCE")),
                        "genomic_locations": {
                            "chr": normalize(docsum.findtext("CHR")),
                            "position": normalize(docsum.findtext("CHRPOS")),
                        },
                        "fxn_class": [],
                        "genes": [],
                        "maf": {},
                        "spdi": normalize(docsum.findtext("SPDI")),
                        "last_updated": normalize(docsum.findtext("UPDATEDATE")),
                    }

                    # Alleles (from DOCSUM SEQ=[C/-] or <ALLELE> tags)
                    alleles_text = docsum.findtext("DOCSUM")
                    if alleles_text and "SEQ=" in alleles_text:
                        seq_part = alleles_text.split("SEQ=")[1].split("|")[0]
                        record["alleles"] = seq_part.strip("[]").split("/")

                    # Functional class
                    fxn = docsum.findtext("FXN_CLASS")
                    if fxn:
                        record["fxn_class"] = [x.strip() for x in fxn.split(",") if x.strip()]

                    # Genes
                    for gene_node in docsum.findall(".//GENES/GENE_E"):
                        record["genes"].append({
                            "symbol": normalize(gene_node.findtext("NAME")),
                            "gene_id": normalize(gene_node.findtext("GENE_ID"))
                        })

                    # Global MAF
                    for maf in docsum.findall(".//GLOBAL_MAFS/MAF"):
                        study = normalize(maf.findtext("STUDY"))
                        freq = normalize(maf.findtext("FREQ"))
                        record["maf"][study] = freq

                    results.append(record)
                except Exception:
                    results.append({
                        "variation_id": normalize(rsid),
                        "accession": f"rs{rsid}",
                        "title": "not specified",
                        "variant_type": "not specified",
                        "alleles": [],
                        "clinical_significance": "not specified",
                        "genomic_locations": {"chr": "not specified", "position": "not specified"},
                        "fxn_class": [],
                        "genes": [],
                        "maf": {},
                        "spdi": "not specified",
                        "last_updated": "not specified",
                    })

            return results

        else:
            raise ValueError(f"Unsupported variant source: {source}")

        return results


    # ---------- Gene Regulation ----------
    def get_regulation(
        self,
        region: str,
        source: str = "encode",
        min_results: int = 3,
        organism: str = "Homo sapiens",
    ) -> List[Dict[str, Any]]:
        """
        Fetch regulatory elements from ENCODE (default) or JASPAR.
        Accepts gene name (converted to genomic region) or explicit region.
        - source: 'encode' (default) or 'jaspar'
        - min_results: number of results per category (default=3, max=20)
        """

        limit = min(20, max(1, min_results))
        results = []

        # ---- ENCODE ----
        if source.lower() == "encode":
            encode_url = f"{self.config.encode_base}/search/"
            categories = {
                "regulatory_elements": "candidate Cis-Regulatory Elements",
                "tf_binding_sites": "footprints",
                "enhancers": "enhancer predictions",
                "promoters": "transcription start sites",
            }

            for label, annot_type in categories.items():
                params = {
                    "type": "Annotation",
                    "annotation_type": annot_type,
                    "organism.scientific_name": organism,
                    "format": "json",
                    "limit": limit,
                }

                try:
                    r_enc = self.session.get(
                        encode_url, params=params, timeout=self.config.timeout
                    )
                    r_enc.raise_for_status()
                except requests.exceptions.HTTPError as e:
                    if r_enc.status_code == 404:
                        # No data for this organism/category → skip gracefully
                        continue
                    else:
                        raise  # re-raise other errors

                enc_data = r_enc.json()
                for hit in enc_data.get("@graph", [])[:limit]:
                    results.append({
                        "source": "ENCODE",
                        "category": label,
                        "id": hit.get("@id", "not specified"),
                        "accession": hit.get("accession", "not specified"),
                        "annotation_type": (
                            hit.get("annotation_type")
                            or hit.get("annotation_subtype")
                            or hit.get("assay_term_name")
                            or annot_type
                        ),
                        "biosample": hit.get("biosample_ontology", {}).get("term_name", "not specified")
                                    if isinstance(hit.get("biosample_ontology"), dict)
                                    else "not specified",
                        "target": hit.get("targets", [{}])[0].get("label", "not specified")
                                    if hit.get("targets") else "not specified",
                        "organism": hit.get("organism", {}).get("scientific_name", "not specified")
                                    if isinstance(hit.get("organism"), dict)
                                    else "not specified",
                        "description": hit.get("description", "not specified"),
                    })

        # ---- JASPAR ----
        elif source.lower() == "jaspar":
            jaspar_api = "https://jaspar.elixir.no/api/v1/matrix"
            params_j = {"search": region, "limit": limit}
            r_j = self.session.get(jaspar_api, params=params_j, timeout=self.config.timeout)
            r_j.raise_for_status()
            jdata = r_j.json()

            for motif in jdata.get("results", [])[:limit]:
                results.append({
                    "source": "JASPAR",
                    "motif_id": motif.get("matrix_id", "not specified"),
                    "name": motif.get("name", "not specified"),
                    "species": motif.get("species", "not specified"),
                    "collection": motif.get("collection_name", "not specified"),
                })

        else:
            raise ValueError(f"Unsupported regulation source: {source}")

        return results

    # ------------- Gene → Disease Associations ----------------
    def get_disease_associations(
        self,
        gene_name: str,
        source: str = "gwas",
        min_results: int = 5,
        organism: str = "Homo sapiens",
    ) -> List[Dict[str, Any]]:
        """
        Fetch gene–disease associations from various sources.
        - gene_name: symbol (e.g., TP53)
        - source: disgenet | gwas | omim | pubmed
        - min_results: number of associations to return (default=5, max=20)
        - organism: organism (default: Homo sapiens)
        """

        results: List[Dict[str, Any]] = []
        limit = min(20, max(1, min_results))

        # ---------- GWAS Catalog ----------
        if source.lower() == "gwas":
            try:
                url = f"https://www.ebi.ac.uk/gwas/rest/api/genes/{gene_name}/associations"
                headers = {"Accept": "application/json"}
                r = self.session.get(url, headers=headers, timeout=self.config.timeout)
                r.raise_for_status()
                data = r.json()

                for assoc in data.get("_embedded", {}).get("associations", [])[:limit]:
                    disease = assoc.get("diseaseTrait", {}).get("trait", "not specified")
                    pub = assoc.get("publicationInfo", {}).get("pubmedId", "not specified")
                    results.append({
                        "source": "GWAS Catalog",
                        "disease": disease,
                        "publications": [pub] if pub != "not specified" else [],
                        "clinical_relevance": assoc.get("studyType", "not specified"),
                    })

            except Exception as e:
                print(f"Error fetching from GWAS Catalog: {e}")

        # ---------- OMIM ----------
        elif source.lower() == "omim":
            try:
                url = "https://api.omim.org/api/entry/search"
                params = {
                    "search": gene_name,
                    "format": "json",
                    "apiKey": self.config.omim_api_key  # must be set in config
                }
                r = self.session.get(url, params=params, timeout=self.config.timeout)
                r.raise_for_status()
                data = r.json()

                entries = (
                    data.get("omim", {})
                        .get("searchResponse", {})
                        .get("entryList", [])
                )

                for entry in entries[:limit]:
                    e = entry.get("entry", {})
                    refs = e.get("referenceList", [])
                    pmids = [ref.get("pubmedID") for ref in refs if "pubmedID" in ref]
                    results.append({
                        "source": "OMIM",
                        "disease": e.get("title", "not specified"),
                        "publications": pmids,
                        "clinical_relevance": e.get("description", "not specified"),
                    })

            except Exception as e:
                print(f"Error fetching from OMIM: {e}")

        # ---------- PubMed ----------
        elif source.lower() == "pubmed":
            try:
                url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                params = {
                    "db": "pubmed",
                    "term": f"{gene_name}[Gene] AND disease",
                    "retmode": "json",
                    "retmax": limit,
                }
                r = self.session.get(url, params=params, timeout=self.config.timeout)
                r.raise_for_status()
                ids = r.json().get("esearchresult", {}).get("idlist", [])

                for pmid in ids:
                    results.append({
                        "source": "PubMed",
                        "disease": "not specified",
                        "publications": [pmid],
                        "clinical_relevance": "Reported association in PubMed",
                    })

            except Exception as e:
                print(f"Error fetching from PubMed: {e}")

        else:
            raise ValueError(f"Unsupported disease association source: {source}")

        return results
