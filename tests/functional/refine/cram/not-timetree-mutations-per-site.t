Setup

  $ source "$TESTDIR"/_setup.sh

Run refine without inferring a time tree.
This is one way to get named internal nodes for downstream analyses and does not require an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations-per-site &> /dev/null

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" \
  >   "$TESTDIR/../data/not_time_tree.nwk" \
  >   tree.nwk \
  >   --significant-digits 2
  {}
  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   "$TESTDIR/../data/mutations_per_site_branch_lengths.json" \
  >   branch_lengths.json \
  >   --significant-digits 0 \
  >   --exclude-paths "root['generated_by']['version']" "root['input_tree']"
  {}
