

from starterator.phamerator_api import fetch_genes_by_pham

print("if you see this i am running")

gaps = fetch_genes_by_pham("1616")

print("number of genes found:", len(gaps))
print()

for gene_name in sorted(gaps):
    print(f"{gene_name}: {gaps[gene_name]}")


"""
from starterator.phams import Pham

pham = Pham("9443")

for gene in pham.genes.values():
    print(gene.full_name)
    
    """