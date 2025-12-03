import dendropy
from dendropy.calculate import treecompare

# Create a shared TaxonNamespace
taxa = dendropy.TaxonNamespace()

# Load both trees using the same TaxonNamespace
t_fast = dendropy.Tree.get(path="fasttree_utf8.nwk", schema="newick", taxon_namespace=taxa)
t_iq = dendropy.Tree.get(path="iqtree_on_prelim-aligned-sequences.treefile", schema="newick", taxon_namespace=taxa)

# Compute unweighted RF distance
rf_distance = treecompare.symmetric_difference(t_fast, t_iq)
print("RF distance (topology only):", rf_distance)