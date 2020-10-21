Integration tests for augur refine.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Try building a time tree.

  $ ${AUGUR} refine \
  >  --tree "refine/tree_raw.nwk" \
  >  --alignment "refine/aligned.fasta" \
  >  --metadata "refine/metadata.tsv" \
  >  --output-tree "$TMP/tree.nwk" \
  >  --output-node-data "$TMP/branch_lengths.json" \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 > /dev/null
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "refine/tree.nwk" "$TMP/tree.nwk" --significant-digits 2
  {}

Build a time tree with mutations as the reported divergence unit.

  $ ${AUGUR} refine \
  >  --tree "refine/tree_raw.nwk" \
  >  --alignment "refine/aligned.fasta" \
  >  --metadata "refine/metadata.tsv" \
  >  --output-tree "$TMP/tree.nwk" \
  >  --output-node-data "$TMP/branch_lengths.json" \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations > /dev/null
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "refine/tree.nwk" "$TMP/tree.nwk" --significant-digits 2
  {}

Run refine without inferring a time tree.
This is one way to get named internal nodes for downstream analyses and does not require an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "refine/tree_raw.nwk" \
  >  --metadata "refine/metadata.tsv" \
  >  --output-tree "$TMP/tree.nwk" \
  >  --output-node-data "$TMP/branch_lengths.json" \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations-per-site > /dev/null
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "refine/not_time_tree.nwk" "$TMP/tree.nwk" --significant-digits 2
  {}
  $ diff -u "refine/mutations_per_site_branch_lengths.json" "$TMP/branch_lengths.json"

Run refine again without a time tree, but request number of mutations per branch as the divergence unit.
This approach only works when we provide an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "refine/tree_raw.nwk" \
  >  --alignment "refine/aligned.fasta" \
  >  --metadata "refine/metadata.tsv" \
  >  --output-tree "$TMP/tree.nwk" \
  >  --output-node-data "$TMP/branch_lengths.json" \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations > /dev/null
  */treetime/aa_models.py:108: VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray (glob)
    [0.800038509951648, 1.20274751778601, 1.55207513886163, 1.46600946033173, 0.830022143283238, 1.5416250309563, 1.53255698189437, 1.41208067821187, 1.47469999960758, 0.351200119909572, 0.570542199221932, 1.21378822764856, 0.609532859331199, 0.692733248746636, 1.40887880416009, 1.02015839286433, 0.807404666228614, 1.268589159299, 0.933095433689795]

Confirm that trees match expected topology and branch lengths, given that the output should not be a time tree.

  $ python3 "$TESTDIR/../../scripts/diff_trees.py" "refine/not_time_tree.nwk" "$TMP/tree.nwk" --significant-digits 2
  {}
  $ diff -u "refine/integer_branch_lengths.json" "$TMP/branch_lengths.json"

Run refine again without a time tree, but try to request number of mutations per branch as the divergence unit.
This approach does not make sense and should not work without an alignment FASTA.

  $ ${AUGUR} refine \
  >  --tree "refine/tree_raw.nwk" \
  >  --metadata "refine/metadata.tsv" \
  >  --output-tree "$TMP/tree.nwk" \
  >  --output-node-data "$TMP/branch_lengths.json" \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 \
  >  --divergence-units mutations > /dev/null
  *ERROR: alignment is required* (glob)
  [1]

  $ popd > /dev/null
