Setup

  $ source "$TESTDIR"/_setup.sh

Confirm that filtering omits strains without metadata or sequences.
The input sequences are missing one strain that is in the metadata.
The metadata are missing one strain that has a sequence.
The list of strains to include has one strain with no metadata/sequence and one strain with information that would have been filtered by country.
The query initially filters 3 strains from Colombia, one of which is added back by the include.

  $ ${AUGUR} filter \
  >  --sequence-index "$TESTDIR/../data/sequence_index.tsv" \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --query "country != 'Colombia'" \
  >  --non-nucleotide \
  >  --exclude-ambiguous-dates-by year \
  >  --include "$TESTDIR/../data/include.txt" \
  >  --output-strains filtered_strains.txt \
  >  --output-log filtered_log.tsv
  4 strains were dropped during filtering
  	1 had no metadata
  	1 had no sequence data
  	3 were filtered out by the query: "country != 'Colombia'"
  \\t1 was added back because it was in .*include\.txt.* (re)
  9 strains passed all filters

  $ head -n 1 filtered_log.tsv; tail -n +2 filtered_log.tsv | sort -k 1,1
  strain	filter	kwargs
  COL/FLR_00008/2015	filter_by_query	"[[""query"", ""country != 'Colombia'""]]"
  COL/FLR_00008/2015\tforce_include_strains\t"[[""include_file"", ""*/data/include.txt""]]" (esc) (glob)
  COL/FLR_00024/2015	filter_by_query	"[[""query"", ""country != 'Colombia'""]]"
  Colombia/2016/ZC204Se	filter_by_query	"[[""query"", ""country != 'Colombia'""]]"
  HND/2016/HU_ME59	filter_by_sequence_index	[]
