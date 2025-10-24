Test `augur filter` with `--metadata-id-columns` and `--include-where`

  $ cat > metadata.tsv <<EOF
  > strain	PPX_accession_version	clade	lineage
  > strain1	PP_000T4L2.1	A	A.1
  > strain2	PP_000T4J6.1	B	B.1
  > strain3	PP_000T4J7.1	B	B.1
  > EOF

  $ cat > sequences.fasta <<EOF
  > >strain1
  > ATCG
  > >strain2
  > ATCG
  > >strain3
  > ATCG
  > EOF

  $ augur filter \
  >   --metadata metadata.tsv \
  >   --metadata-id-columns PPX_accession_version \
  >   --sequences sequences.fasta \
  >   --include-where PPX_accession_version=PP_000T4L2.1 PPX_accession_version=PP_000T4J6.1 \
  >   --output-sequences output-sequences.fasta \
  >   --output-metadata output-metadata.tsv
  4 strains were dropped during filtering
  \t3 had no metadata (esc)
  \t3 had no sequence data (esc)
  \t1 was added back because of 'PPX_accession_version=PP_000T4L2.1' (esc)
  \t1 was added back because of 'PPX_accession_version=PP_000T4J6.1' (esc)
  2 strains passed all filters

  $ cat output-metadata.tsv
  strain	PPX_accession_version	clade	lineage
  strain1	PP_000T4L2.1	A	A.1
  strain2	PP_000T4J6.1	B	B.1

  $ cat output-sequences.fasta
  >strain1
  ATCG
  >strain2
  ATCG

