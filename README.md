# Explore_Biological_Databases

XploreBD is an MCP-compatible extension that provides access to biological data through standardized tools. It integrates APIs such as Ensembl, NCBI, GEO, and ClinVar to enable gene annotation, variant analysis, gene–disease associations, RNA expression queries, and more—making it easy for LLMs to explore biomedical information in a structured and reproducible way.

## Installation

Install the package in development mode:

```bash
pip install -e .
```

Or install from PyPI (when available):

```bash
pip install Xplore-BD
```

## Usage


### Command Line Interface

The package provides a `Xplore-BD` command for usage:

```bash
# Get help
Xplore-BD --help

# Run your main functionality
Xplore-BD <args>
```



### MCP Server

The package provides an MCP server for integration with MCP-compatible clients:

```bash
# Run the MCP server
Xplore-BD-server
```

The MCP server provides the following tools:

- **tool1**: Description of tool1
- **tool2**: Description of tool2
- **tool3**: Description of tool3


### Python API

```python
from Xplore-BD.main import Explore_Biological_DatabasesBridge

# Initialize the bridge
bridge = Explore_Biological_DatabasesBridge()

# Use your functionality
result = bridge.your_method()
```

## Features

- **Feature 1**: Description of feature 1
- **Feature 2**: Description of feature 2
- **Feature 3**: Description of feature 3

- **MCP Integration**: Full Model Context Protocol server implementation


## API Methods

### Core Methods

- `method1()`: Description of method1
- `method2()`: Description of method2
- `method3()`: Description of method3

### Configuration

The package uses a configuration class for settings:

```python
from Xplore-BD.main import Explore_Biological_DatabasesConfig, Explore_Biological_DatabasesBridge

config = Explore_Biological_DatabasesConfig(
    base_url="https://api.example.com",
    api_key="your_api_key",
    timeout=30.0
)

bridge = Explore_Biological_DatabasesBridge(config)
```


## MCP Server Configuration

To use the MCP server with an MCP client, configure it as follows:

```json
{
  "mcpServers": {
    "Xplore-BD": {
      "command": "Xplore-BD-server",
      "env": {}
    }
  }
}
```

The server will automatically handle:
- JSON-RPC communication
- Tool discovery and invocation
- Error handling and reporting


## Development

### Setup Development Environment

```bash
# Install in development mode with dev dependencies
pip install -e .[dev]

# Run tests
pytest

# Format code
black Xplore-BD/

# Type checking
mypy Xplore-BD/
```

### Project Structure

```
Xplore-BD/
├── pyproject.toml      # Package configuration
├── README.md          # This file
├── LICENSE            # MIT License
├── Xplore-BD/         # Main package
│   ├── __init__.py    # Package initialization
│   ├── main.py        # Core functionality

│   └── cli.py         # Command-line interface


│   └── mcp_server.py  # MCP server implementation

└── tests/             # Test files
    ├── __init__.py
    └── test_main.py   # Tests for main functionality
```

## License

MIT License - see LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Run the test suite
6. Submit a pull request

## Support

For issues and questions, please use the GitHub issue tracker. 