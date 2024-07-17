Setup

  $ source "$TESTDIR"/_setup.sh

--output-group-by-sizes does not work without --group-by-weights.

  $ ${AUGUR} filter \
  >   --metadata "$TESTDIR/../data/metadata.tsv" \
  >   --group-by year month \
  >   --subsample-max-sequences 100 \
  >   --output-group-by-sizes target_group_sizes.tsv \
  >   --output-strains strains.txt
  ERROR: --output-group-by-sizes is only available for --group-by-weights. It may be added to other sampling methods in the future - see <https://github.com/nextstrain/augur/issues/new>
  [2]
