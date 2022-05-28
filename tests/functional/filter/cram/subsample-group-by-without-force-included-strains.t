Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Do not consider force-included strains for subsampling.
In this test, we force-include two old strains that are the only representatives of their month/year date group (December 2015).
We don't filter these strains, so they could be considered for subsampling, but Augur removes them from consideration because they have been force-included.

Pandas engine
-------------

  $ cat >$TMP/include_old_strains.txt <<~~
  > PRVABC59
  > COL/FLR_00008/2015
  > ~~

  $ ${AUGUR} filter \
  >   --metadata filter/data/metadata.tsv \
  >   --include $TMP/include_old_strains.txt \
  >   --group-by month year \
  >   --subsample-max-sequences 10 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  $ rm -f $TMP/metadata-filtered.tsv

SQLite engine
-------------

  $ cat >$TMP/include_old_strains.txt <<~~
  > PRVABC59
  > COL/FLR_00008/2015
  > ~~

  $ ${AUGUR} filter --engine sqlite \
  >   --metadata filter/data/metadata.tsv \
  >   --include $TMP/include_old_strains.txt \
  >   --group-by month year \
  >   --subsample-max-sequences 10 \
  >   --output-metadata $TMP/metadata-filtered.tsv > /dev/null
  $ rm -f $TMP/metadata-filtered.tsv
