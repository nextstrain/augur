SETUP

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"

Seqkit errors messages are shown directly.

  $ cat >x.fasta <<~~
  > >seq1
  > ATCG
  > >seq2
  > GCTA
  > >seq3
  > TCGA
  > ~~

  $ cat >y.fasta <<~~
  > invalid fasta file
  > ~~

  $ ${AUGUR} merge \
  >   --sequences x.fasta y.fasta \
  >   --skip-input-sequences-validation \
  >   --output-sequences - > merged.fasta
  Merging sequences and writing to '-'…
  [SeqKit] [ERRO] fastx: invalid FASTA/Q format
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]

Input file doesn't exist

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --output-sequences -
  Validating 'x.fasta'…
  Validating 'z.fasta'…
  ERROR: No such file or directory: 'z.fasta'
  ERROR: Validation failed for 'z.fasta'.
  [2]

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --skip-input-sequences-validation \
  >   --output-sequences -
  Merging sequences and writing to '-'…
  ERROR: No such file or directory: 'z.fasta'
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]
