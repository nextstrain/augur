Integration tests for augur frequencies.

  $ pushd "$TESTDIR" > /dev/null

Generate mutation frequencies grouped by region.

  $ ../../bin/augur frequencies \
  >  --method diffusion \
  >  --metadata frequencies/metadata.tsv \
  >  --alignments frequencies/aligned.fasta \
  >  --gene-names nuc \
  >  --group-by region \
  >  --output "$TMP/mutation_frequencies_by_region.json"

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "frequencies/mutation_frequencies_by_region.json" "$TMP/mutation_frequencies_by_region.json"
  {}

Generate mutation frequencies grouped by country.

  $ augur frequencies \
  >  --method diffusion \
  >  --metadata frequencies/metadata.tsv \
  >  --alignments frequencies/aligned.fasta \
  >  --gene-names nuc \
  >  --group-by country \
  >  --output "$TMP/mutation_frequencies_by_country.json"

  $ python3 "$TESTDIR/../../scripts/diff_jsons.py" "frequencies/mutation_frequencies_by_country.json" "$TMP/mutation_frequencies_by_country.json"
  {}

Try to generate mutation frequencies with an invalid group field.

  $ augur frequencies \
  >  --method diffusion \
  >  --metadata frequencies/metadata.tsv \
  >  --alignments frequencies/aligned.fasta \
  >  --gene-names nuc \
  >  --group-by missing_field \
  >  --output "$TMP/mutation_frequencies.json"
  Error!
