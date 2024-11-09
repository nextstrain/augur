Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with conflicting default arguments.
Expect error message.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "data/aligned.fasta" \
  >  --tree-builder-args="--threads-max 1 --msa data/aligned.fasta" \
  >  --output "tree_raw.nwk"
  ERROR: The following tree builder arguments conflict with hardcoded defaults. Remove these arguments and try again: --threads-max, --msa
  [1]
