Setup

  $ source "$TESTDIR"/_setup.sh

The id column can be used in --query, --exclude-where, and --include-where.

  $ cat > metadata.tsv <<~~
  > strain	accession	clade
  > strain1	acc1	A
  > strain2	acc2	A
  > strain3	acc3	B
  > strain4	acc4	B
  > strain4	acc5	C
  > ~~

  $ ${AUGUR} filter \
  >   --metadata metadata.tsv \
  >   --metadata-id-columns accession \
  >   --query "accession != 'acc5'" \
  >   --exclude-where accession=acc2 clade=B \
  >   --include-where accession=acc3 \
  >   --output-strains filtered.txt
  3 strains were dropped during filtering
  	1 was dropped because of 'accession=acc2'
  	2 were dropped because of 'clade=B'
  	1 was filtered out by the query: "accession != 'acc5'"
  	1 was added back because of 'accession=acc3'
  2 strains passed all filters

  $ sort filtered.txt
  acc1
  acc3
