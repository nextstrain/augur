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
This currently does not work as expected.

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --sequences sequences.fasta \
  >   --output-strains filtered_strains.txt
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  3 strains were dropped during filtering
  \t3 had no sequence data (esc)
  [2]
  $ sort filtered_strains.txt
