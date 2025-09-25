"""
Explore_Biological_Databases - XploreBD is an MCP-compatible extension that provides access to biological data through standardized tools. It integrates APIs such as Ensembl, NCBI, GEO, and ClinVar to enable gene annotation, variant analysis, gene–disease associations, RNA expression queries, and more—making it easy for LLMs to explore biomedical information in a structured and reproducible way.
"""

__version__ = "0.1.0"
__author__ = "Shrihari Kamalan Kumaraguruparan"
__email__ = "kkshrihari@gmail.com"

from .bridge import Bridge, Config

__all__ = ["Bridge", "Config"] 