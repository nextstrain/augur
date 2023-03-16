Integration tests for augur refine.

  $ export AUGUR="${AUGUR:-$TESTDIR/../../bin/augur}"

Try building a time tree.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/refine/tree_raw.nwk" \
  >  --alignment "$TESTDIR/refine/aligned.fasta" \
  >  --metadata "$TESTDIR/refine/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 &> /dev/null

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "$TESTDIR/refine/tree.nwk" tree.nwk --significant-digits 2
  {}

Build a time tree with mutations as the reported divergence unit.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/refine/tree_raw.nwk" \
  >  --alignment "$TESTDIR/refine/aligned.fasta" \
  >  --metadata "$TESTDIR/refine/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations &> /dev/null

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "$TESTDIR/refine/tree.nwk" tree.nwk --significant-digits 2
  {}

Run refine without inferring a time tree.
This is one way to get named internal nodes for downstream analyses and does not require an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/refine/tree_raw.nwk" \
  >  --metadata "$TESTDIR/refine/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations-per-site &> /dev/null

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" \
  >   "$TESTDIR/refine/not_time_tree.nwk" \
  >   tree.nwk \
  >   --significant-digits 2
  {}
  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" \
  >   "$TESTDIR/refine/mutations_per_site_branch_lengths.json" \
  >   branch_lengths.json \
  >   --significant-digits 0 \
  >   --exclude-paths "root['generated_by']['version']" "root['input_tree']"
  {}

Run refine again without a time tree, but request number of mutations per branch as the divergence unit.
This approach only works when we provide an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/refine/tree_raw.nwk" \
  >  --alignment "$TESTDIR/refine/aligned.fasta" \
  >  --metadata "$TESTDIR/refine/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations &> /dev/null

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" \
  >   "$TESTDIR/refine/not_time_tree.nwk" \
  >   tree.nwk \
  >   --significant-digits 2
  {}
  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" \
  >   "$TESTDIR/refine/integer_branch_lengths.json" \
  >   branch_lengths.json \
  >   --significant-digits 0 \
  >   --exclude-paths "root['generated_by']['version']" "root['input_tree']" "root['alignment']"
  {}

Run refine again without a time tree, but try to request number of mutations per branch as the divergence unit.
This approach does not make sense and should not work without an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/refine/tree_raw.nwk" \
  >  --metadata "$TESTDIR/refine/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations > /dev/null
  *ERROR: alignment is required* (glob)
  [1]
