"""Optional external annotation helpers (HGNC and NCBI PubMed)."""
from __future__ import annotations
import time
import urllib.parse
import requests

def get_hgnc_data(gene, cache=None, timeout=5):
    cache = cache if cache is not None else {}
    if gene in cache:
        return cache[gene]
    try:
        response = requests.get(
            f"https://rest.genenames.org/fetch/symbol/{urllib.parse.quote(str(gene))}",
            headers={"Accept": "application/json"}, timeout=timeout
        )
        response.raise_for_status()
        payload = response.json()
        docs = payload.get("response", {}).get("docs", [])
        if docs:
            doc = docs[0]
            cache[gene] = (
                doc.get("symbol", gene),
                doc.get("alias_symbol", []),
                doc.get("prev_symbol", [])
            )
            return cache[gene]
    except (requests.RequestException, ValueError, KeyError, TypeError):
        pass
    cache[gene] = (gene, [], [])
    return cache[gene]

def build_pubmed_query(gene, aliases, previous_symbols):
    terms = sorted(set([str(gene)] + list(aliases) + list(previous_symbols)))
    gene_block = " OR ".join(f'"{term}"' for term in terms)
    return (
        f"({gene_block})[tiab] AND "
        '(azoospermia OR "spermatogenic arrest" OR "meiotic arrest" OR '
        '"Sertoli cell-only" OR hypospermatogenesis OR spermatogenesis OR '
        '"ovarian insufficiency" OR infertil* OR steril*) AND '
        "(mutation OR variant OR deficiency OR knockout OR pathogenic) "
        "NOT (flagell* OR motil*) NOT review[pt]"
    )

def pubmed_score(count):
    if count == 0: return 0
    if count <= 2: return 1
    if count <= 4: return 2
    if count <= 10: return 3
    if count <= 20: return 4
    return 5

def compute_pubscore_strg(symbols, request_delay=0.2, timeout=8, session=None):
    """Query HGNC aliases and PubMed; returns {gene_symbol: score}.

    External service failures are represented by score 0 and do not stop analysis.
    """
    session = session or requests
    hgnc_cache, results = {}, {}
    for symbol in sorted(set(str(s) for s in symbols if s is not None)):
        if not symbol or symbol.lower() == "nan":
            continue
        approved, aliases, previous = get_hgnc_data(symbol, hgnc_cache, timeout=timeout)
        query = build_pubmed_query(approved, aliases, previous)
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        try:
            response = session.get(url, params={"db": "pubmed", "term": query, "retmode": "json"}, timeout=timeout)
            response.raise_for_status()
            count = int(response.json().get("esearchresult", {}).get("count", 0))
            results[symbol] = pubmed_score(count)
        except (requests.RequestException, ValueError, TypeError, KeyError):
            results[symbol] = 0
        if request_delay:
            time.sleep(request_delay)
    return results
