Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	country	date
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > ~~

Create FASTA file for testing.

  $ cat >sequences.fasta <<~~
  > >sequence_A
  > ATCE
  > >sequence_B
  > TCGA
  > >sequence_C
  > TCGA
  > ~~

  $ ${AUGUR} filter \
  >  --sequences sequences.fasta \
  >  --metadata metadata.tsv \
  >  --non-nucleotide \
  >  --output-strains filtered_strains_aa.txt
  1 strain was dropped during filtering
  	1 was dropped because it had non-nucleotide characters (esc)
  2 strains passed all filters
