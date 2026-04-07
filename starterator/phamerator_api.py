import os
import requests
from typing import Any, Dict

#test import

def _auth() -> tuple[str, str]:
    username = os.environ.get("PHAMERATOR_USERNAME")
    password = os.environ.get("PHAMERATOR_PASSWORD")
    if not username or not password:
        raise RuntimeError("Missing PHAMERATOR_USERNAME / PHAMERATOR_PASSWORD env vars")
    return username, password


def _get_json(url: str, timeout: int = 30) -> Any:
    """
    Makes an authenticated GET request and returns parsed JSON.
    No caching. No processing. No storage.
    """
    username, password = _auth()

    print("[phamerator] GET", url, flush=True)
    r = requests.get(url, auth=(username, password), timeout=timeout)
    print("[phamerator] status", r.status_code, flush=True)

    if r.status_code in (401, 403):
        raise RuntimeError("Phamerator auth failed (401/403)")
    if not (200 <= r.status_code < 300):
        raise RuntimeError(f"Phamerator API error {r.status_code}: {(r.text or '')[:200]}")

    return r.json()


def fetch_genes_by_phage(
    phage_name: str,
    dataset: str = "Actino_Draft",
    base_url: str = "https://phamerator.org",
) -> list[dict[str, Any]]:
    """
    Fetch all genes for one phage.

    Endpoint:
      https://phamerator.org/api/<dataset>/genes/<phage_name>

    Returns:
      A list of dictionaries parsed from the API JSON.
    """
    url = f"{base_url.rstrip('/')}/api/{dataset}/genes/{phage_name}"
    data = _get_json(url)

    if not isinstance(data, list):
        raise RuntimeError(f"Expected a list from Phamerator, got {type(data).__name__}")

    return data


def fetch_genes_by_pham(
    pham_no: int | str,
    dataset: str = "Actino_Draft",
    base_url: str = "https://phamerator.org",
) -> Any:
    """
    Fetch data for one pham:

    first fetch gene # for all the genes in the phams

    Then in call the jsons for each of those specific genes in the pham


    """
    url = f"{base_url.rstrip('/')}/api/{dataset}/phamily/{pham_no}"
    return _get_json(url)




def all_gaps_by_phageID(phageID: str) -> Dict[str, int]:
# add function that returns only the gap score from a full phage data
# input is phageID
# returns dictionary, each key is geneID and the value is the gap score

# call fetch_genes_by_phage and then parse the json file

# ***ignore for now

    return None


def all_gaps_by_pham(phamID: str) -> Dict[str, int]:

# returns all gap scores for a pham
# input is phamID
# returns dictionary with key is GeneID and value is gap score

    pham_gaps = {}
    pham_genes = fetch_genes_by_pham(phamID)
    return pham_gaps

test commmit