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
  
  ERROR: Shell exited 255 when running: 
  .*rmdup | (re)
  .*augur write-file - (re)
  
  Command output was:
    \x1b[31m[ERRO]\x1b[0m fastx: invalid FASTA/Q format (esc)
  
  WARNING: Metadata merge was successful but sequence merge has failed, so --output-metadata has no corresponding sequence output.
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]

A missing input file shows a meaningful error during validation…

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --output-sequences -
  Validating 'x.fasta'…
  Validating 'z.fasta'…
  
  ERROR: Shell exited 255 when running: 
  .*fx2tab --name.* (re)
  
  Command output was:
    \x1b[31m[ERRO]\x1b[0m stat z.fasta: no such file or directory (esc)
  
  WARNING: Metadata merge was successful but sequence merge has failed, so --output-metadata has no corresponding sequence output.
  ERROR: Unable to read 'z.fasta'. See error above.
  [2]

… as well as during merging.

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --skip-input-sequences-validation \
  >   --output-sequences -
  Merging sequences and writing to '-'…
  
  ERROR: Shell exited 255 when running: 
  .*rmdup | (re)
  .*augur write-file - (re)
  
  Command output was:
    \x1b[31m[ERRO]\x1b[0m stat z.fasta: no such file or directory (esc)
  
  WARNING: Metadata merge was successful but sequence merge has failed, so --output-metadata has no corresponding sequence output.
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]

At least two sequence inputs are required.

  $ ${AUGUR} merge \
  >   --sequences x.fasta \
  >   --output-sequences -
  WARNING: Metadata merge was successful but sequence merge has failed, so --output-metadata has no corresponding sequence output.
  ERROR: At least two sequence inputs are required for merging.
  [2]
