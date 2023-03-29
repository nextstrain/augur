Setup

  $ source "$TESTDIR"/_setup.sh

Comma-delimited metadata is allowed. However, the output metadata will be tab-delimited.

  $ cat >metadata.txt <<~~
  > strain,column
  > SEQ_1,A
  > SEQ_2,B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt > /dev/null
  $ cat filtered.txt
  strain\tcolumn (esc)
  SEQ_2\tB (esc)

Colon-delimited metadata is allowed. However, the output metadata will be tab-delimited.

  $ cat >metadata.txt <<~~
  > strain:column
  > SEQ_1:A
  > SEQ_2:B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt > /dev/null
  $ cat filtered.txt
  strain\tcolumn (esc)
  SEQ_2\tB (esc)
