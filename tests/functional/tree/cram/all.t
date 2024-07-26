Integration tests for augur tree.

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Try building a tree with IQ-TREE.

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null

Try building a tree with IQ-TREE with more threads (4) than there are input sequences (3).

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 4 > /dev/null

Try building a tree with IQ-TREE using its ModelTest functionality, by supplying a substitution model of "auto".

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --method iqtree \
  >  --substitution-model auto \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null

Build a tree with excluded sites using a compressed input file.

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fasta.xz" \
  >  --exclude-sites "$TESTDIR/../data/excluded_sites.txt" \
  >  --output tree_raw.nwk &> /dev/null

Build a tree, augmenting existing default arguments with custom arguments.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --tree-builder-args="-czb" \
  >  --output tree_raw.nwk > /dev/null

Build a tree, replacing existing default arguments with custom arguments.
Since the following custom arguments are incompatible with the default IQ-TREE arguments, this command will only work with the `--override-default-args` flag.

  $ ${AUGUR} tree \
  >  --method iqtree \
  >  --alignment "$TESTDIR/../data/full_aligned.fasta" \
  >  --tree-builder-args="-czb -bb 1000 -bnni" \
  >  --override-default-args \
  >  --output tree_raw.nwk > /dev/null

Build a tree with an input file that doesn't end in .fasta, and ensure it's not overwritten.

  $ ${AUGUR} tree \
  >  --alignment "$TESTDIR/../data/aligned.fa" \
  >  --method iqtree \
  >  --output tree_raw.nwk \
  >  --nthreads 1 > /dev/null

  $ sha256sum "$TESTDIR/../data/aligned.fa" | awk '{print $1}'
  169a9f5f70b94e26a2c4ab2b3180d4b463112581438515557a9797adc834863d
