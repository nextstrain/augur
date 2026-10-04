Setup

  $ source "$TESTDIR"/_setup.sh

Unaligned sequences to align against the dataset's reference.

  $ python3 -c '
  > from Bio import SeqIO
  > records = list(SeqIO.parse("'"$EBOLA"'/masked.fasta", "fasta"))[:3]
  > for r in records: r.seq = r.seq.replace("-", "")
  > SeqIO.write(records, "sequences.fasta", "fasta")'

The dataset's reference is used as --reference-sequence.

  $ ${AUGUR} align \
  >  --sequences sequences.fasta \
  >  --nextclade-dataset "$EBOLA" \
  >  --output aligned_dataset.fasta > /dev/null 2>&1

  $ ${AUGUR} align \
  >  --sequences sequences.fasta \
  >  --reference-sequence "$EBOLA/reference.fasta" \
  >  --output aligned_explicit.fasta > /dev/null 2>&1

  $ diff aligned_dataset.fasta aligned_explicit.fasta

  $ grep -c ">" aligned_dataset.fasta
  4

--nextclade-dataset and --reference-sequence are mutually exclusive.

  $ ${AUGUR} align \
  >  --sequences sequences.fasta \
  >  --nextclade-dataset "$EBOLA" \
  >  --reference-sequence "$EBOLA/reference.fasta" \
  >  --output aligned.fasta
  ERROR: --nextclade-dataset and --reference-sequence can not be used together.
  [2]
