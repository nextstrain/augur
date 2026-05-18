Setup

  $ source "$TESTDIR"/_setup.sh

Build a tree with conflicting default arguments.
Expect error message.

  $ cp "$TESTDIR/../data/aligned.fasta" .
  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment aligned.fasta \
  >  --tree-builder-args="--threads-max 1 --msa aligned.fasta" \
  >  --output "tree_raw.nwk"
  ERROR: The following tree builder arguments conflict with hardcoded defaults. Remove these arguments and try again: --threads-max, --msa
  [1]
