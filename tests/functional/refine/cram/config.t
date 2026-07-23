Setup

  $ source "$TESTDIR"/_setup.sh

Test that config files work properly with YAML syntax using underscores.

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

  $ ${AUGUR} refine --config config.yaml > /dev/null 2> error.log || cat error.log

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" "$TESTDIR/../data/tree.nwk" tree.nwk --significant-digits 2
  {}

Test that specifying an option in both the config and CLI raises an error.

  $ ${AUGUR} refine --config config.yaml --coalescent opt 2>&1 | grep "error:"
  augur refine: error: The following option(s) were specified both in the config file and on the CLI: coalescent
  [2]
