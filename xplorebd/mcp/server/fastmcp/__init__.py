"""FastMCP - A more ergonomic interface for MCP servers."""
try:
    from importlib.metadata import version
    __version__ = version("mcp")
except Exception:
    __version__ = "0.0.0-local"

from .server import Context, FastMCP
from .utilities.types import Image

__version__ = version("mcp")
__all__ = ["FastMCP", "Context", "Image"]
