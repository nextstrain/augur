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
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  
  Command output was:
    \x1b[31m[ERRO]\x1b[0m fastx: invalid FASTA/Q format (esc)
  
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]

Input file doesn't exist

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --output-sequences -
  Validating 'x.fasta'…
  Validating 'z.fasta'…
  
  ERROR: Shell exited 2 when running: 
  .* (re)
  .* (re)
  
  Command output was:
    ERROR: No such file or directory: 'z.fasta'
    [INFO]\x1b[0m 0 duplicated records removed (esc)
  
  ERROR: Validation failed for 'z.fasta'. See error above.
  [2]

  $ ${AUGUR} merge \
  >   --sequences x.fasta z.fasta \
  >   --skip-input-sequences-validation \
  >   --output-sequences -
  Merging sequences and writing to '-'…
  
  ERROR: Shell exited 2 when running: 
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  .* (re)
  
  Command output was:
    ERROR: No such file or directory: 'z.fasta'
    [INFO]\x1b[0m 0 duplicated records removed (esc)
  
  ERROR: Merging failed, see error(s) above. The command may have already written data to stdout. You may want to clean up any partial outputs.
  [2]
