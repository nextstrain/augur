Setup

  $ source "$TESTDIR"/_setup.sh

Try building a time tree.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 &> /dev/null

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" "$TESTDIR/../data/tree.nwk" tree.nwk --significant-digits 2
  {}
