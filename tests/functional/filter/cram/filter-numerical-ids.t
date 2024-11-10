Setup

  $ source "$TESTDIR"/_setup.sh

Create a pair of files with numerical strain IDs.

  $ cat >metadata.tsv <<~~
  > strain	col1
  > 1	A
  > 2	B
  > 3	C
  > ~~
  $ cat >sequences.fasta <<~~
  > >1
  > AAAA
  > >2
  > AAAA
  > >3
  > AAAA
  > ~~

Test that nothing is filtered out due to missing sequence data.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --output-strains filtered_strains.txt \
  >   2>/dev/null
  $ sort filtered_strains.txt
  1
  2
  3
