Setup

  $ source "$TESTDIR"/_setup.sh

Comma-delimited metadata is allowed by default. However, the output metadata will be tab-delimited.

  $ cat >metadata.txt <<~~
  > strain,column
  > SEQ_1,A
  > SEQ_2,B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt 2>/dev/null
  $ cat filtered.txt
  strain\tcolumn (esc)
  SEQ_2\tB (esc)

Colon-delimited metadata is not allowed by default.

  $ cat >metadata.txt <<~~
  > strain:column
  > SEQ_1:A
  > SEQ_2:B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt > /dev/null
  ERROR: Could not determine the delimiter of 'metadata.txt'. Valid delimiters are: (',', '\t'). This can be changed with --metadata-delimiters.
  [2]

Pass the default valid delimiters explicitly in reverse order.
Note: this shows how to specify a tab character in the list, though it shouldn't be necessary for most users.

  $ cat >metadata.txt <<~~
  > strain:column
  > SEQ_1:A
  > SEQ_2:B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --metadata-delimiters $'\t' ',' \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt > /dev/null
  ERROR: Could not determine the delimiter of 'metadata.txt'. Valid delimiters are: ['\t', ',']. This can be changed with --metadata-delimiters.
  [2]

Allow colon-delimited metadata. However, the output metadata will be tab-delimited.

  $ cat >metadata.txt <<~~
  > strain:column
  > SEQ_1:A
  > SEQ_2:B
  > ~~

  $ ${AUGUR} filter \
  >  --metadata metadata.txt \
  >  --metadata-delimiters ':' \
  >  --exclude-where column=A \
  >  --output-metadata filtered.txt 2>/dev/null
  $ cat filtered.txt
  strain\tcolumn (esc)
  SEQ_2\tB (esc)
