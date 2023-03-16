Setup

  $ source "$TESTDIR"/_setup.sh

Run refine without a time tree, but try to request number of mutations per branch as the divergence unit.
This approach does not make sense and should not work without an alignment FASTA.

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
  >  --divergence-units mutations > /dev/null
  *ERROR: alignment is required* (glob)
  [1]
