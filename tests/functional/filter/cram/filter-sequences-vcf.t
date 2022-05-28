Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Filter TB strains from VCF and save as a list of filtered strains.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --sequences filter/data/tb.vcf.gz \
  >  --metadata filter/data/tb_metadata.tsv \
  >  --min-date 2012 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  $ wc -l "$TMP/filtered_strains.txt"
  \s*3 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --sequences filter/data/tb.vcf.gz \
  >  --metadata filter/data/tb_metadata.tsv \
  >  --min-date 2012 \
  >  --output-strains "$TMP/filtered_strains.txt" > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  $ wc -l "$TMP/filtered_strains.txt"
  \s*3 .* (re)
  $ rm -f "$TMP/filtered_strains.txt"

