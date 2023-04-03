Setup

  $ source "$TESTDIR"/_setup.sh

This file tests the behavior of various mismatched inputs.

Start with a mismatched metadata file and sequence file.

  $ cat >metadata.tsv <<~~
  > strain	col
  > SEQ1	a
  > SEQ2	b
  > SEQ3	c
  > SEQ4	d
  > ~~

  $ cat >partial_sequences.fasta <<~~
  > >SEQ1
  > ATCG
  > >SEQ2
  > ATCG
  > >SEQ5
  > ATCG
  > ~~

Run without a sequence index.

  $ ${AUGUR} filter \
  >  --sequences partial_sequences.fasta \
  >  --metadata metadata.tsv \
  >  --output-strains filtered_strains.txt \
  >  --output-sequences filtered_sequences.fasta \
  >  --output-log filtered_log.txt
  Note: You did not provide a sequence index, so Augur will generate one. You can generate your own index ahead of time with `augur index` and pass it with `augur filter --sequence-index`.
  3 strains were dropped during filtering
  	1 had no metadata
  	2 had no sequence data
  2 strains passed all filters

  $ sort filtered_strains.txt
  SEQ1
  SEQ2

  $ cat filtered_sequences.fasta
  >SEQ1
  ATCG
  >SEQ2
  ATCG

Note that the output says "3 strains were dropped during filtering", but only 2 are shown in the output log.
The one that had no metadata (SEQ5) is not shown here.

  $ tail -n+2 filtered_log.txt | sort
  SEQ3	filter_by_sequence_index	[]
  SEQ4	filter_by_sequence_index	[]

Use a sequence index that has an entry for 1 strain present in both the metadata and FASTA files.

  $ cat >dummy_index_1.tsv <<~~
  > strain	length
  > SEQ1	4
  > ~~

  $ ${AUGUR} filter \
  >  --sequences partial_sequences.fasta \
  >  --sequence-index dummy_index_1.tsv \
  >  --metadata metadata.tsv \
  >  --output-strains filtered_strains.txt \
  >  --output-sequences filtered_sequences.fasta \
  >  --output-log filtered_log.txt
  WARNING: The sequence index is out of sync with the provided sequences. Metadata and strain output may not match sequence output.
  4 strains were dropped during filtering
  	1 had no metadata
  	3 had no sequence data
  1 strains passed all filters

  $ sort filtered_strains.txt
  SEQ1

  $ tail -n+2 filtered_log.txt | sort
  SEQ2	filter_by_sequence_index	[]
  SEQ3	filter_by_sequence_index	[]
  SEQ4	filter_by_sequence_index	[]

  $ cat filtered_sequences.fasta
  >SEQ1
  ATCG

Use a sequence index that has an entry for 1 strain present in both the metadata and FASTA files, and another in just the metadata.

  $ cat >dummy_index_2.tsv <<~~
  > strain	length
  > SEQ1	4
  > SEQ3	4
  > ~~

  $ ${AUGUR} filter \
  >  --sequences partial_sequences.fasta \
  >  --sequence-index dummy_index_2.tsv \
  >  --metadata metadata.tsv \
  >  --output-strains filtered_strains.txt \
  >  --output-sequences filtered_sequences.fasta \
  >  --output-log filtered_log.txt
  WARNING: The sequence index is out of sync with the provided sequences. Metadata and strain output may not match sequence output.
  4 strains were dropped during filtering
  	1 had no metadata
  	2 had no sequence data
  1 strains passed all filters

  $ sort filtered_strains.txt
  SEQ1
  SEQ3

  $ tail -n+2 filtered_log.txt | sort
  SEQ2	filter_by_sequence_index	[]
  SEQ4	filter_by_sequence_index	[]

Note that this doesn't include SEQ3 since it isn't present in the input FASTA file.

  $ cat filtered_sequences.fasta
  >SEQ1
  ATCG
