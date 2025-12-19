import argparse
import dendropy
from dendropy.calculate import treecompare

def get_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("--t1", type=str, required=True, help="Path to first input tree (in Newick format)")
  parser.add_argument("--t2", type=str, required=True, help="Path to first input tree (in Newick format)")
  return parser.parse_args()


def main():
  args = get_args()
    
  # Create a shared TaxonNamespace and load both trees
  taxa = dendropy.TaxonNamespace()
  t_fast = dendropy.Tree.get(path=args.t1, schema="newick", taxon_namespace=taxa)
  t_iq = dendropy.Tree.get(path=args.t2, schema="newick", taxon_namespace=taxa)

  # compute unweighted RF distance
  rf_distance = treecompare.symmetric_difference(t_fast, t_iq)
  print("RF distance (topology only):", rf_distance)

if __name__ == "__main__":
  main()