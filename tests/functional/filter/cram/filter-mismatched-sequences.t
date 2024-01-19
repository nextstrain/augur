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

Run filter. This output makes no sense, but it exists as such because:

- filter_by_sequence_index is not run with --exclude-all.¹
- There is an extra check for valid strains that only applies to this final report, not the actual outputs.

¹ <https://github.com/nextstrain/augur/commit/c256131f93e94eabb91187e33eda2c5736c2bcef>

  $ ${AUGUR} filter \
  >  --sequences sequences.fasta \
  >  --metadata metadata.tsv \
  >  --exclude-all \
  >  --include metadata-ids.txt \
  >  --output-strains filtered_strains.txt
  1 strains were dropped during filtering
  	3 were dropped by `--exclude-all`
  	3 were added back because they were in metadata-ids.txt
  2 strains passed all filters

  $ wc -l filtered_strains.txt
  \s*3 .* (re)
