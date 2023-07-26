Setup

  $ source "$TESTDIR"/_setup.sh

Run export with metadata that contains "accession".

  $ ${AUGUR} export v2 \
  >  --tree "$TESTDIR/../data/tree.nwk" \
  >  --metadata "$TESTDIR/../data/dataset1_metadata_with_strain_and_accession.tsv" \
  >  --node-data "$TESTDIR/../data/div_node-data.json" "$TESTDIR/../data/location_node-data.json" \
  >  --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-accession.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}

Run export with metadata that contains "accession", and use "accession" as the ID column.

  $ ${AUGUR} export v2 \
  >  --tree "$TESTDIR/../data/tree-by-accession.nwk" \
  >  --metadata "$TESTDIR/../data/dataset1_metadata_with_strain_and_accession.tsv" \
  >  --metadata-id-columns accession \
  >  --node-data "$TESTDIR/../data/div_node-data-by-accession.json" "$TESTDIR/../data/location_node-data-by-accession.json" \
  >  --auspice-config "$TESTDIR/../data/auspice_config1.json" \
  >  --maintainers "Nextstrain Team" \
  >  --output dataset.json > /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" "$TESTDIR/../data/dataset-with-accession-by-accession.json" dataset.json \
  >   --exclude-paths "root['meta']['updated']" "root['meta']['maintainers']"
  {}
