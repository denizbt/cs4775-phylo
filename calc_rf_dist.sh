#!/usr/bin/env bash

PY_SCRIPT="rf_distance.py"
TREE_DIR="trees"
OUTFILE="rf_distances.tsv"

printf "group\ttree1\ttree2\trf_distance\n" > "$OUTFILE"

for GROUP in prelim 252 full; do
  case "$GROUP" in
    prelim)
      TREES=(
        "$TREE_DIR/BEAST-prelim-tree.tre"
        "$TREE_DIR/FastTree-prelim-tree.tre"
        "$TREE_DIR/IQ-Tree-prelim-tree.tre"
      )
      ;;
    252)
      TREES=(
        "$TREE_DIR/BEAST-252-tree.tre"
        "$TREE_DIR/FastTree-252-tree.tre"
        "$TREE_DIR/IQ-Tree-252-tree.tre"
        "$TREE_DIR/paper-generated-BEAST-252-tree-ultrafast-bootstrap.tre"
        "$TREE_DIR/paper-generated-BEAST-252-tree.tre"
        "$TREE_DIR/paper-generated-IQ-Tree-252-tree.tre"
      )
      ;;
    full)
      TREES=(
        "$TREE_DIR/BEAST-full-tree.tre"
        "$TREE_DIR/FastTree-full-tree.tre"
        "$TREE_DIR/IQ-Tree-full-tree.tre"
      )
      ;;
  esac

  echo "=== Running group: $GROUP ==="

  for ((i=0; i<${#TREES[@]}; i++)); do
    for ((j=i+1; j<${#TREES[@]}; j++)); do

      t1="${TREES[i]}"
      t2="${TREES[j]}"

      rf=$(python "$PY_SCRIPT" --t1 "$t1" --t2 "$t2" | awk '{print $NF}')

      printf "%s\t%s\t%s\t%s\n" \
        "$GROUP" "$(basename "$t1")" "$(basename "$t2")" "$rf" \
        >> "$OUTFILE"
    done
  done
done