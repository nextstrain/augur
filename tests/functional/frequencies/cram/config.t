Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > method: kde
  > tree: "$TESTDIR/../data/tree.nwk"
  > metadata: "$TESTDIR/../data/metadata.tsv"
  > pivot_interval: 3
  > output: "tip-frequencies.json"
  > ~~

  $ ${AUGUR} frequencies --config config.yaml > /dev/null

  $ diff -u --ignore-matching-lines version "$TESTDIR/../data/zika_tip-frequencies.json" tip-frequencies.json
  $ rm -f tip-frequencies.json
