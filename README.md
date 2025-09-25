# XploreBD – Explore Biological Databases

**XploreBD** is an **MCP-compatible server** that integrates diverse biological data sources into a unified exploration hub.  
Instead of navigating multiple specialized resources (UniProt, GenBank, Ensembl, KEGG, STRING, etc.), XploreBD provides a central interface to query **genes, RNAs, proteins, pathways, and interactions** seamlessly.

By wrapping existing APIs into capability-based tools (e.g., `gene annotate`, `gene variants`, `gene regulation`, `gene diseases`, `literature`), XploreBD simplifies access, ensures consistent input/output formats, and provides a flexible foundation for downstream **bioinformatics, epigenetics, multi-omics research, and machine learning workflows**.

---

## ✨ Features

- 🔎 **Literature search** → PubMed, bioRxiv, medRxiv, ORCID  
- 🧬 **Gene queries**:  
  - Annotation (Ensembl, NCBI)  
  - Expression (GTEx, ENA)  
  - Variants (ClinVar, dbSNP)  
  - Regulation (ENCODE, JASPAR)  
  - Disease associations (GWAS, OMIM, DisGeNET, PubMed)  
- ⚡ **MCP integration** → works as an MCP server for LLMs  
- 🖥 **CLI tool** → rich command-line interface  
- 📦 **Python API** → bridges available for programmatic use  

---

## 📦 Installation

Create and activate an environment (example with Python 3.11):

bash
```
sudo apt update
sudo apt install python3.11 python3.11-venv python3.11-dev

python3.11 -m venv Xplore-env
source Xplore-env/bin/activate
pip install --upgrade pip
```
Install dependencies:
```
pip install -e .
pip install -r requirements.txt
```

Or install from PyPI (when available):
```
pip install xplorebd
```
🚀 Usage
CLI Help
python -m xplorebd --help

📚 Literature

Search publications across PubMed, bioRxiv, medRxiv, ORCID:

python -m xplorebd literature -k "CRISPR gene editing" -c 5

🧬 Gene Annotation
python -m xplorebd gene annotate BRCA1 --source ensembl

📊 Gene Expression
python -m xplorebd gene expression TP53 --organism human

🧬 Gene Variants

Default usage (ClinVar JSON):

python -m xplorebd gene variants TP53


ClinVar as table:

python -m xplorebd gene variants TP53 --table


Explicit ClinVar JSON:

python -m xplorebd gene variants TP53 --source clinvar


Explicit ClinVar table:

python -m xplorebd gene variants TP53 --source clinvar --table


Use dbSNP (JSON):

python -m xplorebd gene variants TP53 --source dbsnp


dbSNP with table:

python -m xplorebd gene variants TP53 --source dbsnp --table


Another gene (ClinVar JSON):

python -m xplorebd gene variants BRCA1


Another gene (dbSNP + table):

python -m xplorebd gene variants BRCA1 --source dbsnp --table

⚖️ Gene Regulation

Input can be a gene name or a genomic region.

Gene name examples:

python -m xplorebd gene regulation TP53 --min 1
python -m xplorebd gene regulation MYC --min 5
python -m xplorebd gene regulation BRCA1 --min 3


Region examples:

python -m xplorebd gene regulation chr17:7661779-7687550 --min 1
python -m xplorebd gene regulation chr8:127735434-127742951 --min 3


Sources:
ENCODE (default):
python -m xplorebd gene regulation TP53 --source encode --min 3


JASPAR (motifs):
python -m xplorebd gene regulation TP53 --source jaspar --min 5


Organism:

python -m xplorebd gene regulation TP53 --organism "Mus musculus" --min 2
python -m xplorebd gene regulation TP53 --organism "Rattus norvegicus" --min 2

Min results (1–20):
python -m xplorebd gene regulation TP53 --min 10

🧾 Gene–Disease Associations
python -m xplorebd gene diseases TP53 --source gwas --min 5
python -m xplorebd gene diseases BRCA1 --source disgenet --min 5
python -m xplorebd gene diseases BRCA1 --source omim

🔌 MCP Server
Run the MCP server for integration with MCP-compatible clients:
python xplorebd/mcp_server.py --serve


manifest.json defines one tool:
literature → search publications across PubMed, bioRxiv, medRxiv, ORCID.
Example MCP client config:

{
  "mcpServers": {
    "xplorebd": {
      "command": "python",
      "args": ["xplorebd/mcp_server.py", "--serve"]
    }
  }
}

🐍 Python API
from xplorebd.bridge import GeneBridge, LiteratureBridge

# Gene annotation
gbridge = GeneBridge()
annotation = gbridge.get_annotation("BRCA1")

# Literature search
lbridge = LiteratureBridge()
papers = lbridge.literature_search(keyword="CRISPR", max_results=3)

⚙️ Configuration
Some APIs require authentication tokens.
Set them as environment variables:
export DISGENET_API_KEY="your_disgenet_key"
export OMIM_API_KEY="your_omim_key"

📂 Project Structure
Xplore-BD/
├── LICENSE              # MIT License
├── README.md            # Documentation
├── manifest.json        # MCP manifest
├── pyproject.toml       # Build system + metadata
├── requirements.txt     # Dependencies
├── xplorebd/            # Main package
│   ├── __init__.py
│   ├── __main__.py
│   ├── bridge.py
│   ├── cli.py
│   ├── gene.py
│   ├── literature.py
│   ├── main.py
│   └── mcp_server.py
└── tests/               # Unit tests

📜 License
This project is licensed under the MIT License – see the LICENSE
 file.

👤 Author
Shrihari Kamalan Kumaraguruparan
📧 kkshrihari@gmail.com

🔗 Project Links
📄 Homepage
📦 Repository
🐛 Issues
