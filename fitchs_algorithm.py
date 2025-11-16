#!/usr/bin/env python3
import argparse
from Bio import Phylo, SeqIO


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str, default="fitch") # require one of Fitch or Sankoff
    parser.add_argument("--newick-tree", type=str, required=True, help="Path to the input tree (in Newick format)")
    parser.add_argument("--fasta-file", type=str, required=True, help="Path to the input fasta files")
    # TODO add more arguments for parameters of the parsimony methods

    return parser.parse_args()


#   TODO Fitch algorithm
def fitch_score(clade, char_map):
    return None

def compute_fitch(tree, sequences):
    # run for each sequence and each char in the sequence
    n = len(next(iter(sequences.values()))) # length of a fasta sequence
    total_score = 0
    for pos in range(n):
        char_map = {a: seq[pos] for a, seq in sequences.items()}
        score = fitch_score(tree.root, char_map)
        # total parsimony score is sum of the scores
        total_score += score
    return total_score


#   TODO Sankoff algorithm (finish Fitch first)
def sankoff(tree, chars, cost):
    pass


def main():
    args = get_args()
    tree = Phylo.read(args.newick_tree, "newick")
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.fasta_file, "fasta")}

    if args.method == "fitch":
          score = compute_fitch(tree, seqs)
          print("Total Fitch score =", score)
    elif args.method == "sankoff":
        sankoff(tree)
    else:
        print(f"Method {args.method} is not supported!")



if __name__ == "__main__":
  main()