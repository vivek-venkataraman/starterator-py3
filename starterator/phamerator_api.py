import os
import requests


def fetch_phage_genes(dataset: str, phage: str, base_url: str = "https://phamerator.org"):
    """
    Calls:
      https://phamerator.org/api/<dataset>/genes/<phage>

    *environment variables contian login information
    """
    username = os.environ.get("PHAMERATOR_USERNAME")
    password = os.environ.get("PHAMERATOR_PASSWORD")
    if not username or not password:
        raise RuntimeError("Missing environment variables")

    url = f"{base_url.rstrip('/')}/api/{dataset}/genes/{phage}"

    print("[phamerator] GET", url, flush=True)
    r = requests.get(url, auth=(username, password), timeout=30)
    print("[phamerator] status", r.status_code, flush=True)

    if r.status_code in (401, 403):
        raise RuntimeError("Phamerator auth failed (401/403)")
    if r.status_code < 200 or r.status_code >= 300:
        raise RuntimeError(f"Phamerator API error {r.status_code}: {(r.text or '')[:200]}")

    return r.json()


def gap_map_from_genes(genes_json):
    """
    Turns the gene list into:
      { "Phage_CDS_1": 17, "Phage_CDS_2": 0, ... }
    """
    gap_map = {}
    for g in genes_json:
        gene_id = g.get("geneID") or g.get("name") or g.get("id")
        gap_map[str(gene_id)] = g.get("gap")
    return gap_map