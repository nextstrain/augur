Setup

  $ source "$TESTDIR"/_setup.sh

Error on duplicates in metadata within same chunk.

  $ cat >metadata-duplicates.tsv <<~~
  > strain	date
  > a	2010-10-10
  > a	2010-10-10
  > b	2010-10-10
  > c	2010-10-10
  > d	2010-10-10
  > ~~
  $ ${AUGUR} filter \
  >   --metadata metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 10 \
  >   --output-metadata metadata-filtered.tsv > /dev/null
  ERROR: The following strains are duplicated in .* (re)
  'a'
  [2]
  $ cat metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

Error on duplicates in metadata in separate chunks.

  $ ${AUGUR} filter \
  >   --metadata metadata-duplicates.tsv \
  >   --group-by year \
  >   --sequences-per-group 2 \
  >   --subsample-seed 0 \
  >   --metadata-chunk-size 1 \
  >   --output-metadata metadata-filtered.tsv > /dev/null
  ERROR: The following strains are duplicated in .* (re)
  'a'
  [2]
  $ cat metadata-filtered.tsv
  cat: .*: No such file or directory (re)
  [1]

Error on duplicates in sequences.

  $ cat >metadata.tsv <<~~
  > strain	date
  > a	2010-10-10
  > b	2010-10-10
  > c	2010-10-10
  > d	2010-10-10
  > ~~

  $ cat >sequences.fasta <<~~
  > >a
  > AAAA
  > >a
  > GGGG
  > >c
  > CCCC
  > >c
  > TTTT
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --output-sequences sequences-filtered.fasta
  ERROR: Sequence ids must be unique.
  
  The following ids were were duplicated in 'sequences.fasta':
  
    a
    c
  
  [2]

Error even if the corresponding output is not used.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --output-strains filtered.txt
  ERROR: Sequence ids must be unique.
  
  The following ids were were duplicated in 'sequences.fasta':
  
    a
    c
  
  [2]

If the error is bypassed with --skip-checks, the output will contain duplicates.

  $ AUGUR_DEBUG=1 ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --skip-checks \
  >   --output-sequences sequences-filtered.fasta
  Skipping first pass of sequence file due to --skip-checks.
  Reading metadata from 'metadata.tsv' and writing to .* (re)
  Writing strains to .* (re)
  Reading sequences from 'sequences.fasta' and writing to 'sequences-filtered.fasta'…
  0 strains were dropped during filtering
  4 strains passed all filters

  $ cat sequences-filtered.fasta
  >a
  AAAA
  >a
  GGGG
  >c
  CCCC
  >c
  TTTT
