Setup

  $ source "$TESTDIR"/_setup.sh

Try building a time tree with a fixed clock rate and clock std dev.

  $ ${AUGUR} refine \
  >  --tree "$TESTDIR/../data/tree_raw.nwk" \
  >  --alignment "$TESTDIR/../data/aligned.fasta" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --output-tree tree.nwk \
  >  --output-node-data branch_lengths.json \
  >  --timetree \
  >  --coalescent opt \
  >  --date-confidence \
  >  --date-inference marginal \
  >  --clock-rate 0.0012 \
  >  --clock-std-dev 0.0002 \
  >  --clock-filter-iqd 4 \
  >  --seed 314159 &> /dev/null

Confirm that JSON output does not include information about the clock rate std dev, since it was provided by the user.

  $ grep -A 4 '\"clock\"' branch_lengths.json
    "clock": {
      "intercept": .*, (re)
      "rate": .*, (re)
      "rtt_Tmrca": .* (re)
    },
