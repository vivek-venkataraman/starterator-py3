from starterator.phamerator_api import fetch_genes_by_phage

genes = fetch_genes_by_phage("Badulia")

print("number of genes:", len(genes))
print("type:", type(genes))

if genes:
    first = genes[0]
    print("first gene keys:", list(first.keys()))
    print("first gene geneID:", first.get("geneID"))
    print("first gene name:", first.get("name"))
    print("first gene gap:", first.get("gap"))