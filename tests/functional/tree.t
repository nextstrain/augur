Integration tests for augur tree.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="${AUGUR:-../../bin/augur}"

Try building a tree with IQ-TREE.

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta \
  >  --method iqtree \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 1 > /dev/null

Try building a tree with IQ-TREE with more threads (4) than there are input sequences (3).

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta \
  >  --method iqtree \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 4 > /dev/null

Try building a tree with IQ-TREE using its ModelTest functionality, by supplying a substitution model of "auto".

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta \
  >  --method iqtree \
  >  --substitution-model auto \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 1 > /dev/null

Build a tree with excluded sites using a compressed input file.

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fasta.xz \
  >  --exclude-sites tree/excluded_sites.txt \
  >  --output "$TMP/tree_raw.nwk" &> /dev/null

Build a tree, augmenting existing default arguments with custom arguments.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment tree/aligned.fasta \
  >  --tree-builder-args="--polytomy" \
  >  --output "$TMP/tree_raw.nwk" > /dev/null

Build a tree, replacing existing default arguments with custom arguments.
Since the following custom arguments are incompatible with the default IQ-TREE arguments, this command will only work with the `--override-default-args` flag.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment tree/full_aligned.fasta \
  >  --tree-builder-args="--polytomy -bb 1000 -bnni --redo" \
  >  --override-default-args \
  >  --output "$TMP/tree_raw.nwk" > /dev/null

Build a tree with conflicting default arguments.
Expect error message.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment tree/aligned.fasta \
  >  --tree-builder-args="--threads-max 1 --msa tree/aligned.fasta" \
  >  --output "$TMP/tree_raw.nwk"
  ERROR: The following tree builder arguments conflict with hardcoded defaults. Remove these arguments and try again: --threads-max, --msa
  [1]

Build a tree with an input file that doesn't end in .fasta, and ensure it's not overwritten.

  $ ${AUGUR} tree \
  >  --alignment tree/aligned.fa \
  >  --method iqtree \
  >  --output "$TMP/tree_raw.nwk" \
  >  --nthreads 1 > /dev/null

  $ sha256sum tree/aligned.fa
  169a9f5f70b94e26a2c4ab2b3180d4b463112581438515557a9797adc834863d  tree/aligned.fa

Clean up tree log files.

  $ rm -f tree/*.log
  $ popd > /dev/null
