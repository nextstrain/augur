Setup

  $ source "$TESTDIR"/_setup.sh

Filter TB strains from VCF and save as a list of filtered strains.

  $ ${AUGUR} filter \
  >  --sequences "$TESTDIR/../data/tb.vcf.gz" \
  >  --metadata "$TESTDIR/../data/tb_metadata.tsv" \
  >  --min-date 2012 \
  >  --output filtered.vcf \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: 'vcftools' is not installed! This is required for VCF data. Please see the augur install instructions to install it.
  [2]
  $ wc -l filtered_strains.txt
  wc: filtered_strains.txt: open: No such file or directory
  [1]

  $ wc -l filtered.vcf
  wc: filtered.vcf: open: No such file or directory
  [1]
