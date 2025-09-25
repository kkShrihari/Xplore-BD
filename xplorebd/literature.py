"""
literature.py — Unified literature search across PubMed, bioRxiv/medRxiv, and ORCID.
"""

import requests
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class Config:
    timeout: float = 30.0
    request_delay: float = 0.2
    ncbi_base: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    europepmc_base: str = "https://www.ebi.ac.uk/europepmc/webservices/rest"
    orcid_base: str = "https://pub.orcid.org/v3.0"


def _clean(s: Optional[str]) -> str:
    """Normalize strings and replace None with 'Not available'."""
    if not s:
        return "Not available"
    return re.sub(r"\s+", " ", s.strip()) or "Not available"


def _filter_ids(ids: Dict[str, Optional[str]]) -> Dict[str, str]:
    """Remove None entries and replace with 'Not available' only if present."""
    return {k: (v if v else "Not available") for k, v in ids.items() if v}


class LiteratureBridge:
    """Bridge for literature search across PubMed, bioRxiv/medRxiv, ORCID."""

    def __init__(self, config: Optional[Config] = None):
        self.config = config or Config()
        self.session = requests.Session()
        self.session.headers.update({"User-Agent": "Explore_BioDB/0.1", "Accept": "application/json"})

    def literature_search(
        self,
        keyword: Optional[str] = None,
        researcher: Optional[str] = None,
        sources: List[str] = ["pubmed", "biorxiv", "medrxiv", "orcid"],
        max_results: int = 20,
    ) -> List[Dict[str, Any]]:
        if not keyword and not researcher:
            raise ValueError("Provide at least one of: keyword or researcher")

        results: List[Dict[str, Any]] = []

        for src in sources:
            if src == "pubmed":
                hits = self._search_pubmed(keyword, researcher, max_results)
            elif src == "biorxiv":
                hits = self._search_europe_pmc(keyword, researcher, max_results, "biorxiv")
            elif src == "medrxiv":
                hits = self._search_europe_pmc(keyword, researcher, max_results, "medrxiv")
            elif src == "orcid" and researcher:
                hits = self._search_orcid(researcher, max_results)
            else:
                hits = []

            for h in hits:
                if len(results) < max_results:
                    results.append(h)
                else:
                    break

            if len(results) >= max_results:
                break

        return results

    # ---------- PubMed ----------
    def _search_pubmed(self, keyword: Optional[str], researcher: Optional[str], max_results: int):
        terms = []
        if keyword:
            terms.append(f"({keyword})")
        if researcher:
            terms.append(f'("{researcher}"[Author])')
        term = " AND ".join(terms)

        url = f"{self.config.ncbi_base}/esearch.fcgi"
        params = {"db": "pubmed", "retmode": "json", "retmax": max_results, "term": term}
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        r.raise_for_status()
        ids = r.json().get("esearchresult", {}).get("idlist", [])

        if not ids:
            return []

        # esummary
        url = f"{self.config.ncbi_base}/esummary.fcgi"
        params = {"db": "pubmed", "retmode": "json", "id": ",".join(ids)}
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        r.raise_for_status()
        summaries = r.json().get("result", {})

        # efetch abstracts
        url = f"{self.config.ncbi_base}/efetch.fcgi"
        params = {"db": "pubmed", "retmode": "xml", "id": ",".join(ids)}
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        r.raise_for_status()
        abstracts = {}
        root = ET.fromstring(r.text)
        for art in root.findall(".//PubmedArticle"):
            pmid_el = art.find(".//ArticleIdList/ArticleId[@IdType='pubmed']")
            pmid = pmid_el.text if pmid_el is not None else None
            abst = " ".join([txt.text or "" for txt in art.findall(".//AbstractText")]).strip()
            if pmid:
                abstracts[pmid] = abst if abst else "Not available"

        out = []
        for pmid in ids:
            doc = summaries.get(pmid, {})
            title = _clean(doc.get("title"))
            authors = [a.get("name") for a in doc.get("authors", []) if a.get("name")] or ["Not available"]
            year = None
            if "pubdate" in doc:
                m = re.search(r"(\d{4})", doc["pubdate"])
                if m:
                    year = int(m.group(1))
            doi = None
            for idobj in doc.get("articleids", []):
                if idobj.get("idtype", "").lower() == "doi":
                    doi = idobj.get("value")
                    break
            out.append({
                "source": "pubmed",
                "title": title,
                "authors": authors,
                "year": year if year else "Not available",
                "venue": doc.get("fulljournalname") or "Not available",
                "citation": f"{', '.join(authors[:3])} ({year}). {title}." if year else "Not available",
                "abstract": abstracts.get(pmid, "Not available"),
                "doi": doi if doi else "Not available",
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "Not available",
                "ids": _filter_ids({"pmid": pmid, "pmcid": None, "orcid": None}),
            })
        return out

    # ---------- Europe PMC ----------
    def _search_europe_pmc(self, keyword: Optional[str], researcher: Optional[str], max_results: int, source: str):
        query_parts = []
        if keyword:
            query_parts.append(f'({keyword})')
        if researcher:
            query_parts.append(f'AUTHOR:"{researcher}"')
        query_parts.append(f'SOURCE:"{source}"')
        query = " AND ".join(query_parts)

        url = f"{self.config.europepmc_base}/search"
        params = {"query": query, "format": "json", "resultType": "core", "pageSize": str(max_results)}
        r = self.session.get(url, params=params, timeout=self.config.timeout)
        r.raise_for_status()
        data = r.json()
        hits = data.get("resultList", {}).get("result", [])

        out = []
        for h in hits:
            authors = [a.get("fullName") for a in h.get("authorList", {}).get("author", []) if a.get("fullName")] if h.get("authorList") else []
            out.append({
                "source": source,
                "title": _clean(h.get("title")),
                "authors": authors if authors else ["Not available"],
                "year": int(h["pubYear"]) if h.get("pubYear") else "Not available",
                "venue": source,
                "citation": f"{h.get('title')} ({h.get('pubYear')})" if h.get("title") else "Not available",
                "abstract": _clean(h.get("abstractText") or h.get("abstract")),
                "doi": h.get("doi") if h.get("doi") else "Not available",
                "url": f"https://doi.org/{h.get('doi')}" if h.get("doi") else (h.get("url") or "Not available"),
                "ids": _filter_ids({"pmid": h.get("pmid"), "pmcid": h.get("pmcid"), "orcid": None}),
            })
        return out

    # ---------- ORCID ----------
    def _search_orcid(self, researcher: str, max_results: int):
        url = f"{self.config.orcid_base}/expanded-search/"
        headers = {"Accept": "application/json"}
        params = {"q": researcher, "rows": str(max_results)}
        r = self.session.get(url, params=params, headers=headers, timeout=self.config.timeout)
        r.raise_for_status()
        data = r.json()
        results = data.get("expanded-result", []) or data.get("result", [])

        out = []
        for rec in results[:max_results]:
            orcid_id = rec.get("orcid-id") or rec.get("orcid")
            out.append({
                "source": "orcid",
                "title": "Not available",
                "authors": [researcher],
                "year": "Not available",
                "venue": "Not available",
                "citation": f"ORCID profile for {researcher} ({orcid_id})" if orcid_id else f"ORCID profile for {researcher}",
                "abstract": "Not available",
                "doi": "Not available",
                "url": f"https://orcid.org/{orcid_id}" if orcid_id else "Not available",
                "ids": _filter_ids({"orcid": orcid_id}),
            })
        return out
