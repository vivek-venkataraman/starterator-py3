import os
import requests
from typing import Any, Dict
from collections import defaultdict
from .database import DB, get_db

def _auth_token() -> str:
    token = os.environ.get("PHAMERATOR_API_TOKEN")
    if not token:
        raise RuntimeError("Missing PHAMERATOR_API_TOKEN env var")
    return token


def _get_json(url: str, timeout: int = 30) -> Any:
    """
    Makes an authenticated GET request and returns parsed JSON.
    No caching. No processing. No storage.
    """
    token = _auth_token()

    print("[phamerator] GET", url, flush=True)
    r = requests.get(url, headers={"Authorization": f"Bearer {token}"}, timeout=timeout)
    print("[phamerator] status", r.status_code, flush=True)

    if r.status_code in (401, 403):
        raise RuntimeError("Phamerator auth failed (401/403)")
    if not (200 <= r.status_code < 300):
        raise RuntimeError(f"Phamerator API error {r.status_code}: {(r.text or '')[:200]}")

    return r.json()


def fetch_genes_by_phage(
    phage_name: str,
    dataset: str = "Actino_Draft",
    base_url: str = "https://classic.phamerator.org",
) -> list[dict[str, Any]]:
    """
    Fetch all genes for one phage.

    Endpoint:
      https://classic.phamerator.org/api/<dataset>/genes/<phage_name>

    Note: https://phamerator.org redirects Actino_Draft requests to
    classic.phamerator.org, and requests/curl both drop the Authorization
    header on a host-changing redirect, so base_url must point at
    classic.phamerator.org directly rather than relying on the redirect.

    Returns:
      A list of dictionaries parsed from the API JSON.
    """
    url = f"{base_url.rstrip('/')}/api/{dataset}/genes/{phage_name}"
    data = _get_json(url)

    if not isinstance(data, list):
        raise RuntimeError(f"Expected a list from Phamerator, got {type(data).__name__}")

    return data



def all_gaps_by_phageID(phageID: str) -> Dict[str, int]:
# add function that returns only the gap score from a full phage data
# input is phageID
# returns dictionary, each key is geneID and the value is the gap score

# call fetch_genes_by_phage and then parse the json file

# ***ignore for now

    return None


def _starterator_db():
    """
    Try to use Starterator's existing/shared DB connection first.
    If that is not available, fall back to making a new DB() connection.
    """
    try:
        return get_db()
    except Exception:
        return DB()


def fetch_genes_by_pham(pham_no, dataset="Actino_Draft", base_url="https://classic.phamerator.org"):
    from .phams import Pham

    pham = Pham(str(pham_no))
    pham_gaps = {}

    genes_by_phage = {}

    for gene in pham.genes.values():
        phage_name = gene.phage_name

        if phage_name not in genes_by_phage:
            genes_by_phage[phage_name] = {}

        gene_num = gene.full_name.split("_")[-1]
        api_gene_id = f"{phage_name}_CDS_{gene_num}".lower()

        genes_by_phage[phage_name][api_gene_id] = gene.full_name

    for phage_name, wanted_genes in genes_by_phage.items():
        genes_json = fetch_genes_by_phage(phage_name, dataset=dataset, base_url=base_url)

        for gene_json in genes_json:
            api_gene_id = gene_json.get("geneID", "").lower()

            if api_gene_id in wanted_genes:
                gene_name = wanted_genes[api_gene_id]
                pham_gaps[gene_name] = gene_json.get("gap")

    return pham_gaps



# create oham object and then do get_gebnes, self.genes property