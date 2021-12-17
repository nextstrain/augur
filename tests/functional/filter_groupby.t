Integration tests for grouping features in augur filter.

  $ pushd "$TESTDIR" > /dev/null
  $ export AUGUR="../../bin/augur"

Try simple grouping.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by country year month \
  >  --subsample-max-sequences 10 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  Sampling at 10 per group.
  2 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  10 strains passed all filters
  $ wc -l "$TMP/filtered_strains.txt"
  \s*10 .* (re)
  $ head -n 2 "$TMP/filtered_metadata.tsv"
  strain\tvirus\taccession\tdate\tregion\tcountry\tdivision\tcity\tdb\tsegment\tauthors\turl\ttitle\tjournal\tpaper_url (esc)
  PRVABC59\tzika\tKU501215\t2015-12-XX\tNorth America\tPuerto Rico\tPuerto Rico\tPuerto Rico\tgenbank\tgenome\tLanciotti et al\thttps://www.ncbi.nlm.nih.gov/nuccore/KU501215\tPhylogeny of Zika Virus in Western Hemisphere, 2015\tEmerging Infect. Dis. 22 (5), 933-935 (2016)\thttps://www.ncbi.nlm.nih.gov/pubmed/27088323 (esc)

Try subsample without any groups.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --subsample-max-sequences 10 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  2 strains were dropped during filtering
  \t2 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  10 strains passed all filters
  $ wc -l "$TMP/filtered_strains.txt"
  \s*10 .* (re)
  $ head -n 2 "$TMP/filtered_metadata.tsv"
  strain\tvirus\taccession\tdate\tregion\tcountry\tdivision\tcity\tdb\tsegment\tauthors\turl\ttitle\tjournal\tpaper_url (esc)
  COL/FLR_00024/2015\tzika\tMF574569\t\tSouth America\tColombia\tColombia\tColombia\tgenbank\tgenome\tPickett et al\thttps://www.ncbi.nlm.nih.gov/nuccore/MF574569\tDirect Submission\tSubmitted (28-JUL-2017) J. Craig Venter Institute, 9704 Medical Center Drive, Rockville, MD 20850, USA\thttps://www.ncbi.nlm.nih.gov/pubmed/ (esc)

Try grouping by an unknown column.
This should warn then continue without grouping.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --exclude-where "region=South America" "region=North America" "region=Southeast Asia" \
  >  --include-where "country=Ecuador" \
  >  --group-by invalid \
  >  --subsample-max-sequences 10 \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  WARNING: The specified group-by categories (['invalid']) were not found. No sequences-per-group sampling will be done.
  10 strains were dropped during filtering
  \t6 of these were dropped because of 'region=South America' (esc)
  \t4 of these were dropped because of 'region=North America' (esc)
  \t1 of these were dropped because of 'region=Southeast Asia' (esc)
  \t1 sequences were added back because of 'country=Ecuador' (esc)
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  2 strains passed all filters
  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*3 .* (re)

Try grouping by an unknown column and a valid column.
This should warn then continue with grouping by valid column only.

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by country invalid \
  >  --subsample-max-sequences 10 \
  >  --subsample-seed 314159 \
  >  --output-metadata "$TMP/filtered_metadata.tsv"
  WARNING: Some of the specified group-by categories couldn't be found: invalid
  Filtering by group may behave differently than expected!
  Sampling at 1 per group.
  3 strains were dropped during filtering
  \t3 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  9 strains passed all filters
  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*10 .* (re)

Try grouping with no probabilistic sampling

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by country year month \
  >  --sequences-per-group 1 \
  >  --subsample-seed 314159 \
  >  --no-probabilistic-sampling \
  >  --output-strains "$TMP/filtered_strains.txt" \
  >  --output-metadata "$TMP/filtered_metadata.tsv" > /dev/null
  $ wc -l "$TMP/filtered_strains.txt"
  \s*9 .* (re)
  $ wc -l "$TMP/filtered_metadata.tsv"
  \s*10 .* (re)

Try grouping with year only

  $ ${AUGUR} filter \
  >  --metadata filter/metadata_ambiguous_months.tsv \
  >  --group-by year \
  >  --sequences-per-group 10 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt"
  0 strains were dropped during filtering
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  10 strains passed all filters

Try grouping with month only

  $ ${AUGUR} filter \
  >  --metadata filter/metadata.tsv \
  >  --group-by month \
  >  --sequences-per-group 10 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt"
  2 strains were dropped during filtering
  \t1 were dropped during grouping due to ambiguous year information (esc)
  \t1 were dropped during grouping due to ambiguous month information (esc)
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  10 strains passed all filters

Try grouping with year, month

  $ ${AUGUR} filter \
  >  --metadata filter/metadata_ambiguous_months.tsv \
  >  --group-by year month \
  >  --sequences-per-group 10 \
  >  --subsample-seed 314159 \
  >  --output-strains "$TMP/filtered_strains.txt"
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  10 strains were dropped during filtering
  \t10 were dropped during grouping due to ambiguous month information (esc)
  \t0 of these were dropped because of subsampling criteria, using seed 314159 (esc)
  [1]
