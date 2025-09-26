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

You can install XploreBD in different ways depending on your workflow.  

### Option 1: Using pip + requirements.txt (recommended for development)

```bash
python3.11 -m venv Xplore-env
source Xplore-env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### Option 2: Using pyproject.toml (PEP 621 build system)

```bash
pip install -e .
```

This uses the dependencies and metadata defined in `pyproject.toml`.  

### Option 3: Using Conda environment (alternative)

```bash
conda create -n xplorebd python=3.11
conda activate xplorebd
pip install -r requirements.txt
```

---

## 🚀 Usage

### CLI Help

```bash
python -m xplorebd --help
```

---

## 📚 Literature

Search publications across **PubMed, bioRxiv, medRxiv, ORCID**:

```bash
python -m xplorebd literature -k "CRISPR gene editing" -c 5
```

---

## 🧬 Gene Annotation

```bash
python -m xplorebd gene annotate BRCA1 --source ensembl
```

---

## 📊 Gene Expression

```bash
python -m xplorebd gene expression TP53 --organism human
```

---

## 🧬 Gene Variants

**Default usage (ClinVar JSON):**
```bash
python -m xplorebd gene variants TP53
```

**ClinVar as table:**
```bash
python -m xplorebd gene variants TP53 --table
```

**Explicit ClinVar JSON:**
```bash
python -m xplorebd gene variants TP53 --source clinvar
```

**Explicit ClinVar table:**
```bash
python -m xplorebd gene variants TP53 --source clinvar --table
```

**Use dbSNP (JSON):**
```bash
python -m xplorebd gene variants TP53 --source dbsnp
```

**dbSNP with table:**
```bash
python -m xplorebd gene variants TP53 --source dbsnp --table
```

**Another gene (ClinVar JSON):**
```bash
python -m xplorebd gene variants BRCA1
```

**Another gene (dbSNP + table):**
```bash
python -m xplorebd gene variants BRCA1 --source dbsnp --table
```

---

## ⚖️ Gene Regulation

Input can be a **gene name** or a **genomic region**.

**Gene name examples:**
```bash
python -m xplorebd gene regulation TP53 --min 1
python -m xplorebd gene regulation MYC --min 5
python -m xplorebd gene regulation BRCA1 --min 3
```

**Region examples:**
```bash
python -m xplorebd gene regulation chr17:7661779-7687550 --min 1
python -m xplorebd gene regulation chr8:127735434-127742951 --min 3
```

**Sources:**
- ENCODE (default):
```bash
python -m xplorebd gene regulation TP53 --source encode --min 3
```

- JASPAR (motifs):
```bash
python -m xplorebd gene regulation TP53 --source jaspar --min 5
```

**Organism:**
```bash
python -m xplorebd gene regulation TP53 --organism "Mus musculus" --min 2
python -m xplorebd gene regulation TP53 --organism "Rattus norvegicus" --min 2
```

**Min results (1–20):**
```bash
python -m xplorebd gene regulation TP53 --min 10
```

---

## 🧾 Gene–Disease Associations

```bash
python -m xplorebd gene diseases TP53 --source gwas --min 5
python -m xplorebd gene diseases BRCA1 --source disgenet --min 5
python -m xplorebd gene diseases BRCA1 --source omim
```

---

## 🔌 MCP Server

Run the MCP server for integration with MCP-compatible clients:

```bash
python xplorebd/mcp_server.py --serve
```

`manifest.json` defines one tool:

- **literature** → search publications across PubMed, bioRxiv, medRxiv, ORCID.

Example MCP client config:

```json
{
  "mcpServers": {
    "xplorebd": {
      "command": "/{path}/Xplore-BD/Xplore-env-py311/bin/python",
      "args": ["xplorebd/mcp_server.py", "--serve"]
    }
  }
}
```

---

## 🐍 Python API

```python
from xplorebd.bridge import GeneBridge, LiteratureBridge

# Gene annotation
gbridge = GeneBridge()
annotation = gbridge.get_annotation("BRCA1")

# Literature search
lbridge = LiteratureBridge()
papers = lbridge.literature_search(keyword="CRISPR", max_results=3)
```

---

## ⚙️ Configuration

Some APIs require authentication tokens.  

Set them as environment variables:

```bash
export DISGENET_API_KEY="your_disgenet_key"
export OMIM_API_KEY="your_omim_key"
```

---

## 📂 Project Structure

```
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
```

---

## 📜 License

This project is licensed under the **MIT License** – see the [LICENSE](LICENSE) file.

---

## 👤 Author

- **Shrihari Kamalan Kumaraguruparan**  
  📧 kkshrihari@gmail.com  

---

## 🔗 Project Links

- 📄 [Homepage](https://github.com/kkShrihari/Xplore-BD)  
- 📦 [Repository](https://github.com/kkShrihari/Xplore-BD)  
- 🐛 [Issues](https://github.com/kkShrihari/Xplore-BD/issues)  
