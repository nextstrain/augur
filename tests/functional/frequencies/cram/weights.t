Setup

  $ source "$TESTDIR"/_setup.sh

Create input files.

  $ cat >metadata.tsv <<~~
  > strain	date	region
  > SEQ1	2021-01-01	A
  > SEQ2	2021-01-02	B
  > SEQ3	2021-01-01	C
  > SEQ4	2021-01-02	C
  > SEQ5	2021-02-02	D
  > ~~

  $ cat >weights.json <<~~
  > { "A": 2, "B": 1, "C": 1, "D": 0 }
  > ~~

  $ cat >tree.nwk <<~~
  > (SEQ2,SEQ3,SEQ4,SEQ5)SEQ1;
  > ~~

Weight by region.

  $ ${AUGUR} frequencies \
  >  --method kde \
  >  --tree tree.nwk \
  >  --metadata metadata.tsv \
  >  --weights weights.json \
  >  --weights-attribute region \
  >  --output tip-frequencies.json >/dev/null

  $ cat tip-frequencies.json
  {
    "SEQ2": {
      "frequencies": [
        0.5,
        0.5
      ]
    },
    "SEQ3": {
      "frequencies": [
        0.0,
        0.0
      ]
    },
    "SEQ4": {
      "frequencies": [
        0.5,
        0.5
      ]
    },
    "SEQ5": {
      "frequencies": [
        0.0,
        0.0
      ]
    },
    "generated_by": {
      "program": "augur",
      "version": ".*" (re)
    },
    "pivots": [
      2021.0041,
      2021.2507
    ]
  } (no-eol)
