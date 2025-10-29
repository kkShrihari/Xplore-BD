# === Inject local vendor (dependencies) so Claude can find them ===
import sys
import os
import asyncio

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# === Import MCP + project bridge ===
try:
    from mcp.server.fastmcp import FastMCPServer as FastMCP
except ImportError:
    # fallback for older SDKs
    from mcp.server.fastmcp import FastMCP

from xplorebd.bridge import Bridge

# Initialize bridge and MCP server
bridge = Bridge()
mcp = FastMCP("XploreBD MCP")

# === Literature Tool ===
@mcp.tool()
async def literature(
    keyword: str = None,
    researcher: str = None,
    max_results: int = 20
):
    """
    Search biomedical literature.

    Sources:
      - PubMed
      - bioRxiv
      - medRxiv
      - ORCID (researcher profiles)

    Args:
      keyword: keyword to search publications
      researcher: ORCID researcher name or ID
      max_results: max number of results (default=20)
    """
    return bridge.literature_search(
        keyword=keyword, researcher=researcher, max_results=max_results
    )

# === Unified Gene Tool ===
@mcp.tool()
async def gene(
    action: str,
    gene: str = None,
    region: str = None,
    source: str = None,
    organism: str = "Homo sapiens",
    min_results: int = 5
):
    """
    Unified Gene Exploration Tool.

    Supported actions:
      - annotation → fetch gene annotation (Ensembl / NCBI)
      - expression → GTEx tissue expression (54 tissues), RNA-seq matrices, ENA fastq links
      - variants → ClinVar (disease variants) or dbSNP (population variants)
      - regulation → enhancers, promoters, TF motifs (ENCODE, JASPAR)
      - diseases → associations with diseases (GWAS, DisGeNET, OMIM, PubMed)

    Args:
      action: one of ["annotation", "expression", "variants", "regulation", "diseases"]
      gene: gene symbol (e.g. "TP53") or Ensembl ID (for annotation/expression/variants/diseases)
      region: genomic region string (e.g. "chr17:7668402-7687550") for regulation
      source: data source override ("ensembl", "ncbi", "clinvar", "dbsnp", "encode", "jaspar", "gwas", "omim", "pubmed")
      organism: default "Homo sapiens"
      min_results: number of results (default=5, max=20)

    Databases used:
      - Ensembl, NCBI, GTEx, ENA, ClinVar, dbSNP, ENCODE, JASPAR, GWAS Catalog, DisGeNET, OMIM, PubMed
    """

    if action == "annotation":
        return bridge.gene_annotation(gene, source or "ensembl")

    elif action == "expression":
        return bridge.gene_expression(gene, organism)

    elif action == "variants":
        return bridge.gene_variants(gene, source or "clinvar", min_results)

    elif action == "regulation":
        return bridge.gene_regulation(region or gene, source or "encode", min_results, organism)

    elif action == "diseases":
        return bridge.gene_diseases(gene, source or "gwas", min_results, organism)

    else:
        return {"error": f"Unsupported action: {action}"}

# === Test Tool ===
@mcp.tool()
async def ping(msg: str = "hello"):
    """Ping tool for connectivity testing."""
    return {"reply": f"pong: {msg}"}

# === Entry point ===
if __name__ == "__main__":
    if "--serve" in sys.argv:
        mcp.run()

