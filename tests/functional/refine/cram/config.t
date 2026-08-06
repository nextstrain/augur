Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > tree: "$TESTDIR/../data/tree_raw.nwk"
  > alignment: "$TESTDIR/../data/aligned.fasta"
  > metadata: "$TESTDIR/../data/metadata.tsv"
  > output_tree: tree.nwk
  > output_node_data: branch_lengths.json
  > timetree: true
  > coalescent: opt
  > date_confidence: true
  > date_inference: marginal
  > clock_filter_iqd: 4
  > seed: 314159
  > ~~

  $ ${AUGUR} refine --config config.yaml > /dev/null

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" "$TESTDIR/../data/tree.nwk" tree.nwk --significant-digits 2
  {}
