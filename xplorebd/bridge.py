"""
Bridge module for Explore_Biological_Databases.

This file acts as a central hub: it re-exports the LiteratureBridge
and GeneBridge (and later OntologyBridge, etc.) so other modules can import
from `bridge` without worrying about file structure.
"""

from .literature import LiteratureBridge, Config
from .gene import GeneBridge

# Alias for backward compatibility / consistency
Bridge = LiteratureBridge

__all__ = [
    "Bridge",
    "Config",
    "LiteratureBridge",
    "GeneBridge",
]
