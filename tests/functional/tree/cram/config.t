Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cp "$TESTDIR/../data/aligned.fasta" .
  $ cat > config.yaml <<~~
  > alignment: "aligned.fasta"
  > method: iqtree
  > output: "tree_raw.nwk"
  > nthreads: 1
  > ~~

  $ ${AUGUR} tree --config config.yaml > /dev/null
