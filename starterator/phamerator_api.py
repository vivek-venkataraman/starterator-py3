import os
import requests

# change to only fetch required genes
def fetch_genes_by_phage(phage_name: str, dataset: str = "Actino_Draft", base_url: str = "https://phamerator.org"):
    """
    https://phamerator.org/api/<dataset>/genes/<phage_name>
    PHAMERATOR_USERNAME / PHAMERATOR_PASSWORD env vars
    """
    username = os.environ.get("PHAMERATOR_USERNAME")
    password = os.environ.get("PHAMERATOR_PASSWORD")
    if not username or not password:
        raise RuntimeError("Missing PHAMERATOR_USERNAME / PHAMERATOR_PASSWORD env vars")

    url = f"{base_url.rstrip('/')}/api/{dataset}/genes/{phage_name}"

    print("[phamerator] GET", url, flush=True)
    r = requests.get(url, auth=(username, password), timeout=30)
    print("[phamerator] status", r.status_code, flush=True)

    if r.status_code in (401, 403):
        raise RuntimeError("Phamerator auth failed (401/403)")
    if not (200 <= r.status_code < 300):
        raise RuntimeError(f"Phamerator API error {r.status_code}: {(r.text or '')[:200]}")

#convert json file to dictionary instead of in json format***

    return r.json()




# processing the json data to return genes from the json: takes only 1 json file at a time for a single gene


def fetch_genes_by_pham
    # endpoint is /phamily/pham #

# return as dictionary