# === Inject local vendor (dependencies) so Claude can find them ===
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
# === Import MCP + project bridge ===
from mcp.server.fastmcp import FastMCP
from xplorebd.bridge import Bridge
import asyncio

# Initialize bridge and MCP server
bridge = Bridge()
mcp = FastMCP("XploreBD Literature MCP")

# === Define MCP tool ===

@mcp.tool()
async def literature(keyword: str = None, researcher: str = None, max_results: int = 20) -> list:
    """
    Search publications across PubMed, bioRxiv/medRxiv, and ORCID.

    Args:
        keyword: Free-text search query (e.g. 'CRISPR gene editing').
        researcher: Researcher full name (e.g. 'Jennifer Doudna').
        max_results: Maximum results per source (default: 20).
    """
    return bridge.literature_search(
        keyword=keyword,
        researcher=researcher,
        max_results=max_results
    )

# === Entry point ===
if __name__ == "__main__":
    if "--serve" in sys.argv:
        mcp.run()
