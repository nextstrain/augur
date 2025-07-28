Integration tests for augur mask.

  $ source "$TESTDIR"/_setup.sh

Mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NTGCTN
  $ rm -f "masked.fasta"

Mask a specific list of sites and also mask one base from the beginning and the end.

  $ ${AUGUR} mask \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --mask-sites 3 4 \
  >  --mask-from-beginning 1 \
  >  --mask-from-end 1 \
  >  --output "masked.fasta"
  Removing masked sites from FASTA file.

  $ cat "masked.fasta"
  >sequence_1
  NTNNTN
  $ rm -f "masked.fasta"
