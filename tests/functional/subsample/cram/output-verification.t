Setup

  $ source "$TESTDIR"/_setup.sh

Test output file format and content verification.

  $ cat >output_test.yaml <<~~
  > samples:
  >   south_america_only:
  >     group_by:
  >     - region
  >     subsample_max_sequences: 2
  >     exclude_where:
  >     - region!=South America
  > ~~

Run subsample to generate outputs for verification.

  $ ${AUGUR} subsample \
  >   --metadata "$TESTDIR"/../../filter/data/metadata.tsv \
  >   --sequences "$TESTDIR"/../../filter/data/sequences.fasta \
  >   --config output_test.yaml \
  >   --output-metadata output_metadata.tsv \
  >   --output-sequences output_sequences.fasta \
  >   --subsample-seed 0
  Validating schema of 'output_test.yaml'...
  Sampling at 2 per group.
  10 strains were dropped during filtering
  	6 were dropped because of 'region!=South America'
  	4 were dropped because of subsampling criteria
  2 strains passed all filters
  11 strains were dropped during filtering
  	1 had no metadata
  	12 were dropped by `--exclude-all`
  \\t2 were added back because they were in .*sample_south_america_only.* (re)
  2 strains passed all filters

Verify output sequence file format.

  $ grep -c '^>' output_sequences.fasta
  2

Verify output metadata file format.

  $ head -n 1 output_metadata.tsv | cut -f1-5
  strain\tvirus\taccession\tdate\tregion (esc)

Verify sequence count matches metadata count.

  $ tail -n +2 output_metadata.tsv | wc -l | tr -d ' '
  2

Verify all sequences in metadata are present in FASTA.

  $ tail -n +2 output_metadata.tsv | cut -f1 | sort > metadata_strains.txt
  $ grep '^>' output_sequences.fasta | sed 's/^>//' | sort > fasta_strains.txt
  $ diff metadata_strains.txt fasta_strains.txt

Verify all sequences are from the correct region.

  $ cut -f5 output_metadata.tsv | tail -n +2 | sort | uniq
  South America
