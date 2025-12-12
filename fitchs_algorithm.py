#!/usr/bin/env python3
import argparse
from Bio import Phylo, SeqIO


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", type=str, default="fitch") # require one of Fitch or Sankoff
    parser.add_argument("--newick-tree", type=str, required=True, help="Path to the input tree (in Newick format)")
    parser.add_argument("--fasta-file", type=str, required=True, help="Path to the input fasta files")
    parser.add_argument("--normalize-labels", action="store_true", help="Try to normalize tree tip labels to match FASTA ids (e.g., strip BEAST descriptions)")
    # TODO add more arguments for parameters of the parsimony methods

    return parser.parse_args()


def fitch_score(clade, char_map):
    #collect observed characters (non-missing) among tips in char_map
    observed = set(v for v in char_map.values() if v not in ("-", "?", "N"))
    if not observed:
        #if everything missing
        return 0

    def helper(node):
        #returns (state_set, score)
        if node.is_terminal():
            ch = char_map.get(node.name, "?")
            if ch in ("-", "?", "N"):
                #treat missing as allowing any observed state
                return set(observed), 0
            return {ch}, 0

        total_score = 0
        child_sets = []
        for child in node.clades:
            s, sc = helper(child)
            child_sets.append(s)
            total_score += sc

        if not child_sets:
            #no children = empty set
            return set(), total_score

        #intersection across children (1+ children
        inter = set.intersection(*child_sets)
        if inter:
            return inter, total_score
        else:
            union = set.union(*child_sets)
            return union, total_score + 1

    _, score = helper(clade)
    return score

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


def normalize_tree_labels(tree, seq_ids):
    """
    1. if full tip name is already in seq_ids keep as-is
    2. use  first whitespace-delimited token
    3. try simple accession-like regex match from start of name
    """
    import re

    access_re = re.compile(r"^([A-Za-z0-9_.-]+)")
    seq_set = set(seq_ids)
    unmatched = []
    changed = 0
    for term in tree.get_terminals():
        name = term.name
        if not name:
            unmatched.append(name)
            continue
        if name in seq_set:
            continue
        #tokenization at whitespace
        token = name.split()[0]
        if token in seq_set:
            term.name = token
            changed += 1
            continue
        #regex match
        m = access_re.match(name)
        if m and m.group(1) in seq_set:
            term.name = m.group(1)
            changed += 1
            continue
        #split on comma for differences like accession / description
        token2 = name.split(",")[0]
        if token2 in seq_set:
            term.name = token2
            changed += 1
            continue
        unmatched.append(name)

    #debugging print
    if changed:
        print(f"Normalized {changed} tip labels to match FASTA IDs")
    if unmatched:
        print(f"Warning: {len(unmatched)} tip labels did not match any FASTA ID; e.g., {unmatched[:5]}")
    return unmatched

def main():
    args = get_args()
    tree = Phylo.read(args.newick_tree, "newick")
    seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.fasta_file, "fasta")}
    if getattr(args, "normalize_labels", False):
        normalize_tree_labels(tree, seqs.keys())

    if args.method == "fitch":
          score = compute_fitch(tree, seqs)
          print("Total Fitch score =", score)
    # elif args.method == "sankoff":
    #     sankoff(tree)
    else:
        print(f"Method {args.method} is not supported!")



if __name__ == "__main__":
  main()