Setup

  $ source "$TESTDIR"/_setup.sh

Calculate diffusion-based tip frequencies from a refined tree with `--regions`.

  $ ${AUGUR} frequencies \
  >  --method diffusion \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --regions "global" "North America" "South America" \
  >  --pivot-interval 3 \
  >  --output tip-frequencies.json > /dev/null

  $ cat tip-frequencies.json
  {
    "BRA/2016/FC_6706": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.2,
        0.2,
        0.2,
        0.2
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "COL/FLR_00008/2015": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.2,
        0.2,
        0.2,
        0.2
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "Colombia/2016/ZC204Se": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.2,
        0.2,
        0.2,
        0.2
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "DOM/2016/BB_0183": {
      "North America": [
        0.25,
        0.25,
        0.25,
        0.25
      ],
      "South America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "EcEs062_16": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.2,
        0.2,
        0.2,
        0.2
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "HND/2016/HU_ME59": {
      "North America": [
        0.25,
        0.25,
        0.25,
        0.25
      ],
      "South America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "PAN/CDC_259359_V1_V3/2015": {
      "North America": [
        0.25,
        0.25,
        0.25,
        0.25
      ],
      "South America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "PRVABC59": {
      "North America": [
        0.25,
        0.25,
        0.25,
        0.25
      ],
      "South America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "VEN/UF_1/2016": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.2,
        0.2,
        0.2,
        0.2
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "ZKC2/2016": {
      "North America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "South America": [
        0.0,
        0.0,
        0.0,
        0.0
      ],
      "global": [
        0.1,
        0.1,
        0.1,
        0.1
      ]
    },
    "counts": {
      "North America": [
        0,
        2,
        2,
        0
      ],
      "South America": [
        0,
        2,
        3,
        0
      ],
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
