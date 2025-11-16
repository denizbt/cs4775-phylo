import sys
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Merge one FASTA file from each subdirectory into a single FASTA."
    )
    parser.add_argument(
        "root",
        help="Root directory containing subdirectories with FASTA files"
    )
    parser.add_argument(
        "-o", "--output",
        default="sequences.fasta",
        help="Output FASTA file name (default: merged.fasta)"
    )
    args = parser.parse_args()

    root = Path(args.root)
    if not root.is_dir():
        sys.exit(f"ERROR: {root} is not a directory")

    fasta_exts = {".fa", ".fasta", ".fna", ".ffn", ".faa"}

    subdirs = sorted([d for d in root.iterdir() if d.is_dir()])
    fasta_files = []

    for d in subdirs:
        # find fasta files inside this subdirectory
        candidates = [
            f for f in d.iterdir()
            if f.is_file() and f.suffix.lower() in fasta_exts
        ]

        if len(candidates) == 0:
            continue
        if len(candidates) > 1:
            print(f"WARNING: {d} has multiple FASTA files; taking all of them", file=sys.stderr)

        fasta_files.extend(sorted(candidates))

    if not fasta_files:
        sys.exit("ERROR: No FASTA files found in any subdirectories")

    with open(args.output, "w") as out:
        for fasta in fasta_files:
            with open(fasta, "r") as f:
                out.write(f.read())

    print(f"Merged {len(fasta_files)} FASTA file(s) --> {args.output}")

if __name__ == "__main__":
    main()