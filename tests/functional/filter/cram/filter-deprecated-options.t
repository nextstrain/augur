Setup

  $ source "$TESTDIR"/_setup.sh

Create files

  $ cat >metadata.tsv <<~~
  > strain	col1
  > SEQ1	A
  > SEQ2	B
  > ~~

  $ cat >sequences.fasta <<~~
  > >SEQ1
  > AAAA
  > >SEQ2
  > AAAA
  > >SEQ3
  > AAAA
  > ~~

--output alias to --output-sequences

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --output filtered.fasta
  WARNING: --output is deprecated. Use --output-sequences instead.
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  1 strain was dropped during filtering
  	1 had no metadata
  2 strains passed all filters

  $ cat filtered.fasta
  >SEQ1
  AAAA
  >SEQ2
  AAAA

-o alias to --output-sequences

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   -o filtered.fasta
  WARNING: -o is deprecated. Use --output-sequences instead.
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  1 strain was dropped during filtering
  	1 had no metadata
  2 strains passed all filters

  $ cat filtered.fasta
  >SEQ1
  AAAA
  >SEQ2
  AAAA
