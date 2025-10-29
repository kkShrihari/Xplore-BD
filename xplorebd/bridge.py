"""
Bridge module for Explore_Biological_Databases.

This file acts as a central hub: it re-exports the LiteratureBridge
and GeneBridge so other modules can import
from `bridge` without worrying about file structure.
"""

from .literature import LiteratureBridge, Config
from .gene import GeneBridge
from .protein import ProteinBridge


class Bridge:
    """Unified Bridge combining Literature and Gene functionality."""

    def __init__(self):
        self._literature = LiteratureBridge()
        self._gene = GeneBridge()

    # --- Literature ---
    def literature_search(self, keyword: str = None, researcher: str = None, max_results: int = 20):
        """Search publications across PubMed, bioRxiv, medRxiv, and ORCID."""
        return self._literature.literature_search(
            keyword=keyword,
            researcher=researcher,
            max_results=max_results,
        )

    # --- Gene ---
    def gene_annotation(self, gene: str, source: str = "ensembl"):
        """Fetch gene annotation from Ensembl or NCBI."""
        return self._gene.get_annotation(gene, source=source)

    def gene_expression(self, gene: str, organism: str = "Homo sapiens"):
        """Query gene expression data (default organism = human)."""
        return self._gene.get_expression(gene, organism=organism)

    def gene_variants(self, gene: str, source: str = "clinvar", min_results: int = 3):
        """Fetch gene variants (ClinVar or dbSNP)."""
        return self._gene.get_variants(gene, source=source, min_results=min_results)

    def gene_regulation(self, region: str, source: str = "encode", min_results: int = 3, organism: str = "Homo sapiens"):
        """Explore gene regulation from ENCODE or JASPAR."""
        return self._gene.get_regulation(region, source=source, min_results=min_results, organism=organism)

    def gene_diseases(self, gene_name: str, source: str = "gwas", min_results: int = 5, organism: str = "Homo sapiens"):
        """Retrieve gene–disease associations from GWAS, DisGeNET, OMIM, or PubMed."""
        return self._gene.get_disease_associations(gene_name, source=source, min_results=min_results, organism=organism)

    # --- Protein ---
    def protein_annotation(self, protein_id: str):
        """Fetch protein annotation (function, domains, sequence features)."""
        return get_protein_annotation(protein_id)

# Export for compatibility
__all__ = [
    "Bridge",
    "Config",
    "LiteratureBridge",
    "GeneBridge",
    "ProteinBridge",
]
