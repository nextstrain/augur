Setup

  $ source "$TESTDIR"/_setup.sh

'KX369547.1' is removed with --clock-filter-iqd 2.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --timetree \
  >  --clock-filter-iqd 2 \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  pruning leaf  KX369547.1

  $ grep -q -F 'KX369547.1' tree.nwk && echo 'Present' || echo 'Pruned'
  Pruned

Use --keep-ids to force-include it.

  $ cat > include.txt <<~~
  > KX369547.1  # Keep me
  > ~~

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --timetree \
  >  --clock-filter-iqd 2 \
  >  --keep-ids include.txt \
  >  --seed 314159 2>&1 | grep "pruning leaf" || echo "Nothing pruned"
  Nothing pruned

  $ grep -q -F 'KX369547.1' tree.nwk && echo 'Present' || echo 'Pruned'
  Present
