"""
protein.py — Unified protein annotation (UniProt + NCBI + Human Protein Atlas)
"""

import re
import requests
from typing import Dict, Any, Optional


class ProteinBridge:
    """Bridge for fetching protein annotation (UniProt first, fallback to NCBI, plus optional HPA expression)."""

    def __init__(self, user_agent: str = "Explore_BioDB/0.1", timeout: int = 30):
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": user_agent,
            "Accept": "application/json"
        })
        self.timeout = timeout

    # ---------------- INTERNAL UTILITIES ---------------- #

    def _clean_dict(self, d):
        """Recursively remove empty or null values."""
        if isinstance(d, dict):
            return {k: self._clean_dict(v) for k, v in d.items() if v not in ("", None, [], {})}
        elif isinstance(d, list):
            return [self._clean_dict(x) for x in d if x not in ("", None, [], {})]
        else:
            return d

    def _get_length(self, location: dict) -> int:
        """Calculate residue span length (if possible)."""
        try:
            start = location.get("start", {}).get("value")
            end = location.get("end", {}).get("value")
            if start and end:
                return int(end) - int(start)
        except Exception:
            pass
        return 0

    def _resolve_to_uniprot_id(self, query: str) -> Optional[str]:
        """
        Resolve a protein/gene name (like TP53 or BRCA1) to a UniProt accession (e.g. P04637).
        Defaults to Homo sapiens (human) first, but falls back to other species if needed.
        """
        if re.match(r"^[A-Z0-9]{6,10}$", query):
            return query

        search_url = "https://rest.uniprot.org/uniprotkb/search"

        try:
            # 1️⃣ Try reviewed HUMAN first
            params = {
                "query": f"(gene_exact:{query} OR protein_name:{query}) AND reviewed:true AND organism_id:9606",
                "fields": "accession",
                "format": "json",
                "size": 1
            }
            r = self.session.get(search_url, params=params, timeout=self.timeout)
            if r.ok:
                results = r.json().get("results", [])
                if results:
                    return results[0].get("primaryAccession")

            # 2️⃣ Try reviewed (any organism)
            params["query"] = f"(gene_exact:{query} OR protein_name:{query}) AND reviewed:true"
            r = self.session.get(search_url, params=params, timeout=self.timeout)
            if r.ok:
                results = r.json().get("results", [])
                if results:
                    return results[0].get("primaryAccession")

            # 3️⃣ Try unreviewed (TrEMBL)
            params["query"] = f"gene_exact:{query} OR protein_name:{query}"
            r = self.session.get(search_url, params=params, timeout=self.timeout)
            if r.ok:
                results = r.json().get("results", [])
                if results:
                    return results[0].get("primaryAccession")
        except Exception:
            pass

        return None

    # ---------------- HUMAN PROTEIN ATLAS ---------------- #

    # def _get_hpa_expression(self, protein_name: str) -> Optional[Dict[str, Any]]:
    #    """
    #    Fetch expression data (RNA, tissue, cell line, and images)
    #    from the Human Protein Atlas (HPA) for a given protein/gene name.
    #    Parses structured JSON embedded in the main HPA HTML page.
    #    """
    #    import json
    #    import re
    #    try:
    #        # Step 1: resolve Ensembl ID from Ensembl REST API
    #        ens_url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{protein_name}?content-type=application/json"
    #        r = self.session.get(ens_url, timeout=self.timeout)
    #        if not r.ok:
    #            return None
    #        ensembl_hits = r.json()
    #        ensembl_id = None
    #        for h in ensembl_hits:
    #            if h.get("type") == "gene":
    #                ensembl_id = h.get("id")
    #                break
    #        if not ensembl_id:
    #            return None

    #        # Step 2: fetch the HPA HTML page
    #        page_url = f"https://www.proteinatlas.org/{ensembl_id}-{protein_name}"
    #        r = self.session.get(page_url, timeout=self.timeout)
    #        if not r.ok:
    #            return {"gene": protein_name, "ensembl_id": ensembl_id, "source_url": page_url}

    #        html = r.text

    #        # Step 3: extract the JSON block inside <script type="application/ld+json">
    #        match = re.search(r'<script[^>]+type="application/ld\+json"[^>]*>(.*?)</script>', html, re.DOTALL)
    #        if not match:
    #            return {"gene": protein_name, "ensembl_id": ensembl_id, "source_url": page_url}

    #        json_text = match.group(1).strip()
    #        try:
    #            data = json.loads(json_text)
    #        except Exception:
    #            # occasionally, multiple JSONs are concatenated; take the first valid one
    #            json_text = re.split(r'\n\s*\n', json_text)[0]
    #            data = json.loads(json_text)

    #        # Step 4: extract expression-related details
    #        expr = {
    #            "gene": protein_name,
    #            "ensembl_id": ensembl_id,
    #            "summary": data.get("description") or f"{protein_name} expression summary.",
    #            "source_url": page_url,
    #        }

    #        # extract tissue expression (if present)
    #        expr_data = data.get("mainEntity", {}).get("expression", [])
    #        if expr_data:
    #            expr["tissue_expression"] = [
    #                {
    #                    "tissue": e.get("tissue"),
    #                    "level": e.get("level"),
    #                    "reliability": e.get("reliability")
    #                }
    #                for e in expr_data if isinstance(e, dict)
    #            ]

    #        # extract image links
    #        images = []
    #        for img in data.get("image", []):
    #            if isinstance(img, dict):
    #                images.append({
    #                    "caption": img.get("caption"),
    #                    "url": img.get("contentUrl")
    #                })
    #        if images:
    #            expr["images"] = images

    #        return self._clean_dict(expr)

    #    except Exception as e:
    #         print(f"[DEBUG] Failed to fetch HPA data: {e}")
    #         return None


    def get_protein_structure(self, protein_input: str) -> dict:
        """
        Fetch 3D (PDB), homology (SWISS-MODEL), and AlphaFold2 predicted structures.
        Works with protein name, gene symbol, or UniProt ID.
        """
        import requests, re, urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

        s = requests.Session()
        s.headers.update({"User-Agent": "XploreBioDB/1.0", "Accept": "application/json"})

        # --- Detect sequence ---
        def _is_sequence(seq: str) -> bool:
            seq = seq.strip().upper().replace("\n", "").replace(" ", "")
            return len(seq) >= 30 and all(c in "ACDEFGHIKLMNPQRSTVWY" for c in seq)

        # --- Resolve UniProt ID ---
        def _resolve_to_uniprot(query: str):
            query = query.strip()
            if re.match(r"^[A-Z0-9]{6,10}$", query):
                return query, "Unknown organism"

            url = "https://rest.uniprot.org/uniprotkb/search"
            variants = [
                f"(gene_exact:{query} OR protein_name:{query}) AND reviewed:true AND organism_id:9606",
                f"(gene_exact:{query.lower()} OR protein_name:{query.lower()}) AND reviewed:true AND organism_id:9606",
                f"(gene:{query} OR protein_name:{query}) AND reviewed:true",
                query,
            ]
            for v in variants:
                try:
                    r = s.get(url, params={
                        "query": v,
                        "fields": "accession,organism_name",
                        "format": "json",
                        "size": 2,
                    }, timeout=15, verify=False)
                    if r.ok and r.json().get("results"):
                        res = r.json()["results"][0]
                        return (
                            res.get("primaryAccession"),
                            res.get("organism", {}).get("scientificName", "Unknown organism"),
                        )
                except Exception:
                    continue
            return None, None

        # --- Handle sequence input ---
        if _is_sequence(protein_input):
            seq = protein_input.strip().replace("\n", "").replace(" ", "")
            return {
                "input_type": "sequence",
                "uniprot_id": None,
                "organism": None,
                "source": ["AlphaFold DB"],
                "structures": {
                    "AlphaFold_model": {
                        "status": "predicted_from_sequence",
                        "sequence_length": len(seq),
                        "browser_url": "https://alphafold.ebi.ac.uk/",
                        "note": "Predicted directly from sequence (no UniProt ID required)."
                    }
                }
            }

        # --- Resolve UniProt ID ---
        uid, organism = _resolve_to_uniprot(protein_input)
        if not uid:
            raise ValueError(f"Could not resolve '{protein_input}' to a UniProt ID.")

        out = {
            "input_type": "name_or_id",
            "uniprot_id": uid,
            "organism": organism,
            "source": ["RCSB PDB", "AlphaFold DB", "SWISS-MODEL"],
            "structures": {},
        }

        # --- AlphaFold2 model ---
        try:
            r = s.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uid}", timeout=20, verify=False)
            if r.ok and r.text.strip():
                data = r.json()
                if isinstance(data, list) and data:
                    d = data[0]
                    #print(f"[DEBUG] AlphaFold API response for {uid}: {list(d.keys())}")
                    conf_score = d.get("confidenceScore") or d.get("globalMetricValue")
                    #print(f"[DEBUG] AlphaFold confidence for {uid}: {conf_score}")
                    out["structures"]["AlphaFold_model"] = {
                        "entry_id": d.get("entryId", f"AF-{uid}-F1"),
                        "description": d.get("uniprotDescription"),
                        "organism": d.get("organismScientificName"),
                        "confidence_score": conf_score,
                        "pdb_url": d.get("pdbUrl"),
                        "cif_url": d.get("cifUrl"),
                        "pae_image": d.get("paeImageUrl") or d.get("paeDocUrl"),
                        "browser_url": f"https://alphafold.ebi.ac.uk/entry/{uid}",
                        "method": "AlphaFold2 predicted model"
                    }
                else:
                    out["structures"]["AlphaFold_model"] = {
                        "browser_url": f"https://alphafold.ebi.ac.uk/entry/{uid}",
                        "method": "AlphaFold2 model (link only)"
                    }
            else:
                out["structures"]["AlphaFold_model"] = {
                    "browser_url": f"https://alphafold.ebi.ac.uk/entry/{uid}",
                    "method": "AlphaFold2 model (link only)"
                }
        except Exception as e:
            print(f"[DEBUG] AlphaFold fetch failed: {e}")
            out["structures"]["AlphaFold_model"] = {
                "browser_url": f"https://alphafold.ebi.ac.uk/entry/{uid}",
                "method": "AlphaFold2 model (link only)"
            }

        # --- Experimental structures (RCSB PDB) ---
        try:
            query = {
                "query": {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": uid,
                    },
                },
                "return_type": "entry",
            }
            r = s.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query, timeout=20, verify=False)
            if r.ok and r.text.strip():
                entries = []
                for res in r.json().get("result_set", [])[:10]:
                    pid = res["identifier"]
                    meta = s.get(f"https://data.rcsb.org/rest/v1/core/entry/{pid}", timeout=10, verify=False)
                    if not meta.ok:
                        continue
                    m = meta.json()
                    expt = m.get("exptl", [{}])[0]
                    entries.append({
                        "pdb_id": pid,
                        "method": expt.get("method", "Unknown"),
                        "resolution": (
                            m.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
                            or expt.get("resolution") or "Not available"
                        ),
                        "release_date": m.get("rcsb_accession_info", {}).get("initial_release_date", "Unknown"),
                        "rcsb_url": f"https://www.rcsb.org/structure/{pid}",
                    })
                out["structures"]["3D_structure"] = {"entries": entries}
        except Exception as e:
            print(f"[DEBUG] PDB fetch failed: {e}")

        # --- Homology models (SWISS-MODEL Repository) ---
        try:
            sm_url = f"https://swissmodel.expasy.org/repository/uniprot/{uid}.json"
            r = s.get(sm_url, timeout=20, verify=False)
            if r.ok and r.text.strip():
                data = r.json()
                structures = data.get("result", {}).get("structures", [])
                models = []
                for model in structures[:3]:  # limit to top 3
                    template = model.get("template")
                    models.append({
                        "template": template,
                        "coverage": model.get("coverage"),
                        "method": "Homology modeling (SWISS-MODEL)",
                        "coordinates": model.get("coordinates"),
                        "model_page": f"https://swissmodel.expasy.org/repository/uniprot/{uid}?template={template}" if template else None,
                        "created_date": model.get("created_date")
                    })
                out["structures"]["homology_model"] = {
                    "repository_url": f"https://swissmodel.expasy.org/repository/uniprot/{uid}",
                    "entries": models
                }
            else:
                out["structures"]["homology_model"] = {"entries": []}
        except Exception as e:
            print(f"[DEBUG] SWISS-MODEL fetch failed: {e}")
            out["structures"]["homology_model"] = {"entries": []}

        return out
