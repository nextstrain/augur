Setup

  $ source "$TESTDIR"/_setup.sh

Create metadata TSV file for testing.

  $ cat >metadata.tsv <<~~
  > strain	country	date
  > sequence_A	USA	2020-10-01
  > sequence_B	USA	2020-10-02
  > sequence_C	USA	2020-10-03
  > ~~

Create FASTA file for testing.

  $ cat >sequences.fasta <<~~
  > >sequence_A
  > ATCG
  > >sequence_B
  > TCGA
  > ~~

Create a strains file to with all strains from metadata.

  $ cat metadata.tsv | cut -f 1 | tail -n +2 > metadata-ids.txt

Run filter. sequence_C is still output even though it is not in sequences
because --include takes precedence.

  $ ${AUGUR} filter \
  >  --sequences sequences.fasta \
  >  --metadata metadata.tsv \
  >  --exclude-all \
  >  --include metadata-ids.txt \
  >  --output-strains filtered_strains.txt
  0 strains were dropped during filtering
  	3 were dropped by `--exclude-all`
  	3 were added back because they were in metadata-ids.txt
  3 strains passed all filters

  $ wc -l filtered_strains.txt
  \s*3 .* (re)
