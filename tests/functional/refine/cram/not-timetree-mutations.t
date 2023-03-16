Setup

  $ source "$TESTDIR"/_setup.sh

Run refine without a time tree, but request number of mutations per branch as the divergence unit.
This approach only works when we provide an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations &> /dev/null

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" \
  >   "$TESTDIR/../data/not_time_tree.nwk" \
  >   tree.nwk \
  >   --significant-digits 2
  {}
  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/integer_branch_lengths.json" \
  >   branch_lengths.json \
  >   --significant-digits 0 \
  >   --exclude-paths "root['generated_by']['version']" "root['input_tree']" "root['alignment']"
  {}
