Test `augur filter` with `--metadata-id-columns`, `--exclude-where`, and `--include-where`

  $ cat > metadata.tsv <<EOF
  > strain	accession	clade
  > strain1	acc1	A
  > strain2	acc2	A
  > strain3	acc3	B
  > strain4	acc4	B
  > EOF

  $ cat > sequences.fasta <<EOF
  > >acc1
  > ATCG
  > >acc2
  > ATCG
  > >acc3
  > ATCG
  > >acc4
  > ATCG
  > EOF

  $ augur filter \
  >   --metadata metadata.tsv \
  >   --metadata-id-columns accession \
  >   --sequences sequences.fasta \
  >   --exclude-where accession=acc2 clade=B \
  >   --include-where accession=acc3 \
  >   --output-sequences output-sequences.fasta \
  >   --output-metadata output-metadata.tsv
  2 strains were dropped during filtering
  \t1 was dropped because of 'accession=acc2' (esc)
  \t2 were dropped because of 'clade=B' (esc)
  \t1 was added back because of 'accession=acc3' (esc)
  2 strains passed all filters

  $ cat output-metadata.tsv
  strain	accession	clade
  strain1	acc1	A
  strain3	acc3	B

  $ cat output-sequences.fasta
  >acc1
  ATCG
  >acc3
  ATCG

