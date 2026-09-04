Setup

  $ source "$TESTDIR"/_setup.sh

Calculate diffusion-based tip frequencies from a refined tree with
`--minimal-clade-size-to-estimate`

  $ ${AUGUR} frequencies \
  >  --method diffusion \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --pivot-interval 3 \
  >  --minimal-clade-size-to-estimate 5 \
  >  --output tip-frequencies.json > /dev/null

  $ cat tip-frequencies.json
  {
    "BRA/2016/FC_6706": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "COL/FLR_00008/2015": {
      "global": [
        0.108525,
        0.107021,
        0.092939,
        0.089031
      ]
    },
    "Colombia/2016/ZC204Se": {
      "global": [
        0.108525,
        0.107021,
        0.092939,
        0.089031
      ]
    },
    "DOM/2016/BB_0183": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "EcEs062_16": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "HND/2016/HU_ME59": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "PAN/CDC_259359_V1_V3/2015": {
      "global": [
        0.108525,
        0.107021,
        0.092939,
        0.089031
      ]
    },
    "PRVABC59": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "VEN/UF_1/2016": {
      "global": [
        0.108525,
        0.107021,
        0.092939,
        0.089031
      ]
    },
    "ZKC2/2016": {
      "global": [
        0.094317,
        0.09532,
        0.104707,
        0.107313
      ]
    },
    "counts": {
      "global": [
        0,
        5,
        5,
        0
      ]
    },
    "generated_by": {
      "program": "augur",
      "version": ".*" (re)
    },
    "pivots": [
      2015.7521,
      2016.0041,
      2016.2527,
      2016.5014
    ]
  } (no-eol)
