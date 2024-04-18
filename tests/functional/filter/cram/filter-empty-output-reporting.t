Setup

  $ source "$TESTDIR"/_setup.sh

Filter using `--exclude-all` to easily get an empty result.

Test the default behavior for empty results is an error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude-all \
  >  --output-strains filtered_strains.txt > /dev/null
  ERROR: All samples have been dropped! Check filter rules and metadata file format.
  [2]
  $ wc -l filtered_strains.txt
  wc: filtered_strains.txt: open: No such file or directory
  [1]

Repeat with the --empty-output-reporting=warn option.
This whould output a warning message but no error.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --exclude-all \
  >  --output-strains filtered_strains.txt \
  >  --empty-output-reporting warn > /dev/null
  WARNING: All samples have been dropped! Check filter rules and metadata file format.
  $ wc -l filtered_strains.txt
  \s*0 .* (re)

Ignore empty results with the --empty-output-reporting=silent option.
Make sure all 3 output types are empty, except the metadata output should still include the header.
This should not output any messages to stderr.

  $ ${AUGUR} filter \
  >  --metadata "$TESTDIR/../data/metadata.tsv" \
  >  --sequences "$TESTDIR/../data/sequences.fasta" \
  >  --exclude-all \
  >  --output-sequences filtered_seqs.fasta \
  >  --output-metadata filtered_metadata.tsv \
  >  --output-strains filtered_strains.txt \
  >  --empty-output-reporting silent > /dev/null
  $ wc -l filtered_seqs.fasta
  \s*0 .* (re)
  $ diff <(head -n 1 filtered_metadata.tsv) <(head -n 1 "$TESTDIR/../data/metadata.tsv")
  $ wc -l filtered_strains.txt
  \s*0 .* (re)
