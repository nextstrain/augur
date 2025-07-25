Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Try masking a VCF without any specified mask.

  $ ${AUGUR} mask --sequences "$TESTDIR/../data/variants.vcf.gz"
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]

Try masking sequences without any specified mask.

  $ ${AUGUR} mask --sequences "$TESTDIR/../data/sequences.fasta"
  No masking sites provided. Must include one of --mask, --mask-from-beginning, --mask-from-end, --mask-invalid, or --mask-sites
  [1]
