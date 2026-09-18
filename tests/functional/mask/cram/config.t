Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > sequences: "$TESTDIR/../data/sequences.fasta"
  > mask_from_beginning: 1
  > mask_from_end: 1
  > mask_sites: [3, 4]
  > output: "masked.fasta"
  > ~~

  $ ${AUGUR} mask --config config.yaml
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NTNNTN
