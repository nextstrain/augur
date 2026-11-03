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
  Validating schema of 'config.yaml'...

Confirm that TreeTime trees match expected topology and branch lengths.

  $ python3 "$TESTDIR/../../../../scripts/diff_trees.py" "$TESTDIR/../data/tree.nwk" tree.nwk --significant-digits 2
  {}

Paths in the config file that can't be resolved show a clean error.

  $ cat > config-missing-path.yaml <<~~
  > tree: missing.nwk
  > ~~

  $ ${AUGUR} refine --config config-missing-path.yaml
  Validating schema of 'config-missing-path.yaml'...
  ERROR: Path 'missing.nwk' not resolvable from any of the following paths:
  
    .
    .* (re)
  [2]
