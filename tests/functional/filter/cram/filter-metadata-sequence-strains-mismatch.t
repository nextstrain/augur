Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Confirm that filtering omits strains without metadata or sequences.
The input sequences are missing one strain that is in the metadata.
The metadata are missing one strain that has a sequence.
The list of strains to include has one strain with no metadata/sequence and one strain with information that would have been filtered by country.
The query initially filters 3 strains from Colombia, one of which is added back by the include.

Pandas engine
-------------

  $ ${AUGUR} filter \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --query "country != 'Colombia'" \
  >  --non-nucleotide \
  >  --exclude-ambiguous-dates-by year \
  >  --include filter/data/include.txt \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-log "$TMP/filtered_log.tsv"
  4 strains were dropped during filtering
  \t1 had no metadata (esc)
  \t1 had no sequence data (esc)
  \t3 of these were filtered out by the query: "country != 'Colombia'" (esc)
  \t1 strains were added back because they were in filter/data/include.txt (esc)
  9 strains passed all filters

  $ diff -u <(sort -k 1,1 filter/data/filtered_log.tsv) <(sort -k 1,1 "$TMP/filtered_log.tsv")
  $ rm -f "$TMP/filtered_strains.txt"

SQLite engine
-------------

  $ ${AUGUR} filter --engine sqlite \
  >  --sequence-index filter/data/sequence_index.tsv \
  >  --metadata filter/data/metadata.tsv \
  >  --query "country != 'Colombia'" \
  >  --non-nucleotide \
  >  --exclude-ambiguous-dates-by year \
  >  --include filter/data/include.txt \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-log "$TMP/filtered_log.tsv"
  4 strains were dropped during filtering
  \t1 had no metadata (esc)
  \t1 had no sequence data (esc)
  \t3 of these were filtered out by the query: "country != 'Colombia'" (esc)
  \t1 strains were added back because they were in filter/data/include.txt (esc)
  9 strains passed all filters

  $ diff -u <(sort -k 1,1 filter/data/filtered_log.tsv) <(sort -k 1,1 "$TMP/filtered_log.tsv")
  $ rm -f "$TMP/filtered_strains.txt"
