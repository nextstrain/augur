Setup

  $ source "$TESTDIR"/_setup.sh

Do not consider force-included strains for subsampling.
In this test, we force-include two old strains that are the only representatives of their month/year date group (December 2015).
We don't filter these strains, so they could be considered for subsampling, but Augur removes them from consideration because they have been force-included.

  $ cat >include_old_strains.txt <<~~
  > PRVABC59
  > COL/FLR_00008/2015
  > ~~

  $ ${AUGUR} filter \
  >   --metadata "$TESTDIR/../data/metadata.tsv" \
  >   --include include_old_strains.txt \
  >   --group-by month year \
  >   --subsample-max-sequences 10 \
  >   --output-metadata metadata-filtered.tsv 2>/dev/null
