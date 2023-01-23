Setup

  $ source "$TESTDIR"/_setup.sh

Filter TB strains from VCF and save as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/tb.vcf.gz" \
  >  --metadata "$TESTDIR/../data/tb_metadata.tsv" \
  >  --min-date 2012 \
  >  --output-strains filtered_strains.txt > /dev/null
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  $ wc -l filtered_strains.txt
  \s*3 .* (re)
