Setup

  $ source "$TESTDIR"/_setup.sh

Filter TB strains from VCF and save as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/tb.vcf.gz" \
  >  --metadata "$TESTDIR/../data/tb_metadata.tsv" \
  >  --min-date 2012 \
  >  --output-sequences filtered.vcf \
  >  --output-strains filtered_strains.txt > /dev/null
  162 strains were dropped during filtering
  	155 had no sequence data
  	7 were dropped because they were earlier than 2012.0 or missing a date
  3 strains passed all filters
  $ wc -l filtered_strains.txt
  \s*3 .* (re)

  $ wc -l filtered.vcf
  \s*2314 .* (re)
