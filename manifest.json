"""
Command-line interface for Explore_Biological_Databases.
"""

import argparse
import sys
import json
from typing import Optional
from .protein import ProteinBridge
from .bridge import LiteratureBridge, GeneBridge, Config


def create_parser() -> argparse.ArgumentParser:
   parser = argparse.ArgumentParser(
      description="XploreBD - Access to biological databases (Literature, Genes, etc.)",
      formatter_class=argparse.RawDescriptionHelpFormatter,
   )

   subparsers = parser.add_subparsers(dest="command", help="Available commands")

   # ---------- Literature search ----------
   lit_parser = subparsers.add_parser("literature", help="Search publications")
   lit_parser.add_argument("--keyword", "-k", help="Free text search query")
   lit_parser.add_argument(
      "--researcher", "-r",
      nargs="+",   # allows multiple words without quotes
      help="Researcher full name (e.g., Jane Doe)"
   )
   lit_parser.add_argument(
      "--sources", "--source",
      dest="sources",
      type=str.lower,
      help="Comma-separated sources (pubmed,biorxiv,medrxiv,orcid)",
      default="pubmed,biorxiv,medrxiv,orcid",
   )
   lit_parser.add_argument(
      "-c", "--count",
      type=int,
      default=3,
      help="Maximum number of results to fetch (default: 3)"
   )
   lit_parser.add_argument("-o", "--output", help="Write JSON results to file")

   # ---------- Gene queries ----------
   gene_parser = subparsers.add_parser("gene", help="Query gene-related information")
   gene_subparsers = gene_parser.add_subparsers(dest="gene_cmd", help="Gene sub-commands")

   # Gene annotation
   ann_parser = gene_subparsers.add_parser(
      "annotate",
      aliases=["gene_annotate", "GA", "geneann"],
      help="Get gene annotation"
   )
   ann_parser.add_argument("gene_id", help="Gene ID or symbol (e.g., BRCA1, ENSG...)")
   ann_parser.add_argument(
      "--source",
      choices=["ensembl", "ncbi"],
      type=str.lower,
      default="ensembl",
      help="Annotation source (default: ensembl)"
   )

   # Gene expression
   expr_parser = gene_subparsers.add_parser("expression", help="Get gene expression profiles")
   expr_parser.add_argument("gene_name", help="Gene name (e.g., BRCA1)")
   expr_parser.add_argument("--organism", default="human", help="Organism (default: human)")

   # Gene variants
   var_parser = gene_subparsers.add_parser("variants", help="Get gene variants")
   var_parser.add_argument("gene_name", help="Gene name (e.g., TP53)")
   var_parser.add_argument(
      "--source",
      choices=["clinvar", "dbsnp"],
      default="clinvar",
      help="Variant source (default: clinvar)"
   )
   var_parser.add_argument(
      "--min",
      type=int,
      default=3,
      help="Number of variants to return (default: 3, max: 20)"
   )
   var_parser.add_argument("--table", action="store_true", help="Pretty-print as table instead of JSON")
   
    # gene regulation
   reg_parser = gene_subparsers.add_parser(
        "regulation",
        aliases=["R", "reg", "gene_regulation", "genereg"],
        help="Fetch gene regulatory elements (ENCODE default, or JASPAR)"
    )
   reg_parser.add_argument(
        "region",
        help="Gene name (e.g., TP53) or region (e.g., chr17:43044295-43170245)"
    )
   reg_parser.add_argument(
        "--min",
        type=int,
        default=3,
        help="Number of regulatory elements to return (default=3, max=20)"
    )
   reg_parser.add_argument(
        "--source",
        choices=["encode", "jaspar"],
        default="encode",
        help="Source database (default: encode)"
    )
   reg_parser.add_argument(
    "--organism",
    type=str,
    default="Homo sapiens",
    help="Organism to query (default: Homo sapiens)"
    )



    # gene diseases
   dis_parser = gene_subparsers.add_parser(
    "diseases",
    aliases=["disease", "assoc"],
    help="Fetch disease associations for a gene (OMIM, GWAS Catalog, DisGeNET)"
    )
   dis_parser.add_argument(
    "gene_name",
    help="Gene symbol (e.g., TP53)"
    )
   dis_parser.add_argument(
    "--min",
    type=int,
    default=3,
    help="Number of associations to return (default=3, max=20)"
    )
   dis_parser.add_argument(
    "--source",
    choices=["omim", "gwas", "disgenet"],
    default="disgenet",
    help="Source database (default: disgenet)"
    )
   dis_parser.add_argument(
    "--organism",
    type=str,
    default="Homo sapiens",
    help="Organism to query (default: Homo sapiens)"
    )


   # ---------- Protein Annotation ----------
   p_annot_parser = subparsers.add_parser(
       "p_annotation",
       help="Fetch protein functional annotation from UniProt (default) or NCBI"
   )
   p_annot_parser.add_argument(
       "protein_id",
       help="Protein accession or name (e.g., P04637, TP53, NP_000537)"
   )
   p_annot_parser.add_argument(
       "--source",
       choices=["uniprot", "ncbi"],
       default="uniprot",
       help="Annotation source (default: uniprot)"
   )
   p_annot_parser.add_argument(
       "--organism",
       type=str,
       default="Homo sapiens",
       help="Organism to query (default: Homo sapiens)"
   )
   p_annot_parser.add_argument(
       "--reviewed",
       choices=["yes", "no"],
       default="no",
       help="Restrict to reviewed entries if source=uniprot (default: no)"
   )

   # ---------- Protein Expression (Human Protein Atlas only) ----------
   pexpr_parser = subparsers.add_parser(
       "p_expression",
       help="Fetch protein expression data from the Human Protein Atlas (human only)"
   )
   pexpr_parser.add_argument(
       "protein_id",
       help="Protein or gene name (e.g., TP53, BRCA1)"
   )

   # ---------- Protein Structure (PDB + AlphaFold) ----------
   pstruct_parser = subparsers.add_parser(
       "p_structure",
       aliases=["protein_structure"],
       help="Fetch 3D, homology, and AlphaFold models for a protein name, UniProt ID, or raw sequence"
   )
   pstruct_parser.add_argument(
       "protein_input",
       help="Protein name, UniProt ID, or raw amino acid sequence (≥30 aa)"
   )
   pstruct_parser.add_argument(
       "--source",
       choices=["both"],
       default="both",
       help="Select structure source:  both (default: both)"
   )
   pstruct_parser.add_argument(
       "--organism",
       type=str,
       default="Homo sapiens",
       help="Organism to query (default: Homo sapiens)"
   )

    # ---------- Protein–Disease Associations (UniProt + DisGeNET + ClinVar + OMIM) ----------
   pdisease_parser = subparsers.add_parser(
        "p_disease",
        aliases=["protein_disease"],
        help="Fetch disease associations for a given protein or gene from UniProt, DisGeNET, ClinVar, and OMIM."
   )

   pdisease_parser.add_argument(
        "protein_input",
        help="Protein name, gene symbol, or UniProt ID (e.g., TP53, BRCA1, or P04637)."
   )

   pdisease_parser.add_argument(
        "--organism",
        type=str,
        default="human",
        help="Organism to query (default: human; options: human, mouse, rat, zebrafish)."
   )

   return parser


def main(args: Optional[list] = None) -> int:
   parser = create_parser()
   parsed_args = parser.parse_args(args)

   if not parsed_args.command:
      parser.print_help()
      return 1

   try:
      # ---------- Literature ----------
      if parsed_args.command == "literature":
         bridge = LiteratureBridge()
         srcs = [s.strip() for s in parsed_args.sources.split(",") if s.strip()]
         researcher_name = " ".join(parsed_args.researcher) if parsed_args.researcher else None

         results = bridge.literature_search(
            keyword=parsed_args.keyword,
            researcher=researcher_name,
            sources=srcs,
            max_results=parsed_args.count,
         )

         if parsed_args.output:
            with open(parsed_args.output, "w", encoding="utf-8") as f:
               json.dump(results, f, indent=2, ensure_ascii=False)
            print(f"Wrote {len(results)} results to {parsed_args.output}")
         else:
            print(json.dumps(results, indent=2, ensure_ascii=False))

      # ---------- Gene ----------
      elif parsed_args.command == "gene":
         gbridge = GeneBridge()

         if parsed_args.gene_cmd in ["annotate", "gene_annotate", "GA", "geneann"]:
            res = gbridge.get_annotation(parsed_args.gene_id, source=parsed_args.source)
            print(json.dumps(res, indent=2, ensure_ascii=False))

         elif parsed_args.gene_cmd == "expression":
            res = gbridge.get_expression(parsed_args.gene_name, organism=parsed_args.organism)
            print(json.dumps(res, indent=2, ensure_ascii=False))

         elif parsed_args.gene_cmd == "variants":
            # cap at 20
            min_results = min(parsed_args.min, 20)
            res = gbridge.get_variants(
               parsed_args.gene_name,
               source=parsed_args.source,
               min_results=min_results
            )

            if parsed_args.table and isinstance(res, list):
               headers = [
                  "variation_id",
                  "accession",
                  "title",
                  "variant_type",
                  "clinical_significance",
                  "review_status",
                  "last_evaluated",
                  "allele_origin",
                  "cytogenetic_location",
                  "protein_change",
                  "molecular_consequence",
                  "conditions",
                  "submissions",
               ]
               print("\t".join(headers))

               def clean(val):
                  if not val or str(val).lower() in ("unknown", "null", "n/a"):
                     return "not specified"
                  return str(val)

               for r in res:
                  conditions = ";".join(r.get("conditions", []) or []) or "not specified"
                  mc = ";".join(r.get("molecular_consequence", []) or []) or "not specified"
                  subs = ";".join(r.get("submissions", []) or []) or "not specified"

                  row = [
                     clean(r.get("variation_id")),
                     clean(r.get("accession")),
                     clean(r.get("title")),
                     clean(r.get("variant_type")),
                     clean(r.get("clinical_significance")),
                     clean(r.get("review_status")),
                     clean(r.get("last_evaluated")),
                     clean(r.get("allele_origin")),
                     clean(r.get("cytogenetic_location")),
                     clean(r.get("protein_change")),
                     mc,
                     conditions,
                     subs,
                  ]
                  print("\t".join(row))
            else:
               print(json.dumps(res, indent=2, ensure_ascii=False))

         elif parsed_args.gene_cmd in ["regulation", "R", "reg", "gene_regulation", "genereg"]:
            res = gbridge.get_regulation(
            parsed_args.region,
            min_results=parsed_args.min,
            organism=parsed_args.organism
            )
            print(json.dumps(res, indent=2, ensure_ascii=False))


         elif parsed_args.gene_cmd in ["diseases", "disease", "assoc"]:
            res = gbridge.get_disease_associations(
            parsed_args.gene_name,
            source=parsed_args.source,
            min_results=parsed_args.min,
            organism=parsed_args.organism
            )
            print(json.dumps(res, indent=2, ensure_ascii=False))

         else:
            print(f"Unknown command: {parsed_args.command}")
            return 1

      # ---------- Protein Annotation ----------
      elif parsed_args.command == "p_annotation":
         pbridge = ProteinBridge()
         res = pbridge.get_protein_annotation(
            protein_id=parsed_args.protein_id
         )
         print(json.dumps(res, indent=2, ensure_ascii=False))
         return 0

      # ---------- Protein Expression ----------
      elif parsed_args.command == "p_expression":
         pbridge = ProteinBridge()
         res = pbridge.get_expression_only(
            protein_name=parsed_args.protein_id
         )
         print(json.dumps(res, indent=2, ensure_ascii=False))
         return 0

      # ---------- Protein Structure ----------
      elif parsed_args.command in ["p_structure", "protein_structure"]:
         pbridge = ProteinBridge()
         user_input = parsed_args.protein_input.strip()

         try:
            res = pbridge.get_protein_structure(user_input)
            res["organism"] = parsed_args.organism

            # --- Filter by selected source ---
            src = parsed_args.source.lower()

            if src == "alphafold":
                # Show only AlphaFold2 model
                res["structures"] = {
                    "AlphaFold_model": res["structures"].get("AlphaFold_model")
                }

            elif src == "pdb":
                # Show only experimental and homology structures
                res["structures"] = {
                    "3D_structure": res["structures"].get("3D_structure"),
                    "homology_model": res["structures"].get("homology_model")
                }

            elif src == "both":
                # Ensure AlphaFold is included if available
                alpha = res["structures"].get("AlphaFold_model")
                pdb = res["structures"].get("3D_structure")
                homo = res["structures"].get("homology_model")
                res["structures"] = {}
                if alpha:
                    res["structures"]["AlphaFold_model"] = alpha
                if pdb:
                    res["structures"]["3D_structure"] = pdb
                if homo:
                    res["structures"]["homology_model"] = homo

            print(json.dumps(res, indent=2, ensure_ascii=False))
            return 0

         except Exception as e:
            print(f"Error fetching protein structure: {e}", file=sys.stderr)
            return 1

      # ---------- Protein–Disease Associations ----------
      elif parsed_args.command in ["p_disease", "protein_disease"]:
         pbridge = ProteinBridge()
         user_input = parsed_args.protein_input.strip()

         try:
            # --- Fetch disease associations from all databases automatically ---
            res = pbridge.get_protein_diseases(
                protein_input=user_input,
                organism=parsed_args.organism.lower()
            )

            # --- Print formatted JSON output ---
            print(json.dumps(res, indent=2, ensure_ascii=False))
            return 0

         except Exception as e:
            print(f"Error fetching protein disease associations: {e}", file=sys.stderr)
            return 1

   except Exception as e:
      print(f"Error: {e}", file=sys.stderr)
      return 1


if __name__ == "__main__":
   sys.exit(main())
