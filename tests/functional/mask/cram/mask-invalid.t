Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Mask invalid nucleotides

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/invalidnucleotide.fasta" \
  >  --mask-invalid \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  ATCGNNNN
  $ rm -f "masked.fasta"
