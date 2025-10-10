Setup

  $ source "$TESTDIR"/_setup.sh

Try building a time tree with mismatched sequence and metadata IDs.
This should produce a warning because the default 'strain' column has wrong IDs.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata_mismatch.tsv" \
  >  --output-tree tree_mismatch.nwk \
  >  --output-node-data branch_lengths_mismatch.json \
  >  --timetree \
  >  --seed 314159 2>&1 | grep -A 5 "WARNING"

  WARNING: For 11 sequence IDs, only 0 corresponding metadata rows could be matched.
    Metadata has 12 entries in total.
    Metadata is using 'strain' as the ID column.
    You may need to explicitly set the metadata ID column using --metadata-id-columns.
    By default, the columns ('strain', 'name') are tried in order.


Now try with the correct ID column specified.
This should work without the ID mismatch warning.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata_mismatch.tsv" \
  >  --metadata-id-columns name \
  >  --output-tree tree_correct.nwk \
  >  --output-node-data branch_lengths_correct.json \
  >  --timetree \
  >  --seed 314159 2>&1 | grep -c "corresponding metadata rows"
  0
