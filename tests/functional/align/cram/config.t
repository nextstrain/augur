Setup

  $ source "$TESTDIR"/_setup.sh

Run with a config file.

  $ cat > config.yaml <<~~
  > sequences:
  >   - "$TESTDIR/../data/sequences.fasta"
  > output: "aligned.fasta"
  > nthreads: 1
  > ~~

  $ ${AUGUR} align --config config.yaml > /dev/null

  $ cat "aligned.fasta"
  >with_gaps
  ---ATATA---
  >no_gaps
  GGGATATAGGG
  >some_other_seq
  --GATCTAGGG
