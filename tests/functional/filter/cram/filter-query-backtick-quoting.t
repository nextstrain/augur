Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata file for testing.

  $ cat >metadata.tsv <<~~
  > strain	region name
  > SEQ_1	A
  > SEQ_2	A
  > SEQ_3	B
  > SEQ_4	
  > ~~

The 'region name' column is query-able by backtick quoting.

  $ ${AUGUR} filter \
  >  --metadata metadata.tsv \
  >  --query '(`region name` == "A")' \
  >  --output-strains filtered_strains.txt > /dev/null
  WARNING: Could not infer columns from the pandas query. Reading all metadata columns,
  which may impact execution time. If the query is valid, please open a new issue:
  
      <https://github.com/nextstrain/augur/issues/new/choose>
  
  and add the query to the description:
  
      (`region name` == "A")

  $ sort filtered_strains.txt
  SEQ_1
  SEQ_2
