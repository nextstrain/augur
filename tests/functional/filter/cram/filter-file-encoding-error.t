Setup

  $ source "$TESTDIR"/_setup.sh

Create a metadata file that contains a non-ASCII character.

  $ cat >metadata.tsv <<~~
  > strain	col1
  > SEQ_1	ã
  > SEQ_2	b
  > SEQ_3	c
  > ~~

Encode it as WINDOWS-1252.

  $ iconv -f UTF-8 -t WINDOWS-1252 metadata.tsv > metadata-windows-1252.tsv

The UTF-8 encoded file can be used without issues.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --output-strains filtered_strains.txt
  0 strains were dropped during filtering
  3 strains passed all filters

An error is shown when using the WINDOWS-1252 encoded file.

  $ ${AUGUR} filter \
  >  --metadata metadata-windows-1252.tsv \
  >  --output-strains filtered_strains.txt
  ERROR: File 'metadata-windows-1252.tsv' contains b'\xe3' which is not valid in the expected 'utf-8' encoding.
  Try re-saving the file using the 'utf-8' encoding.
  [2]
