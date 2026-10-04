Setup

  $ source "$TESTDIR"/_setup.sh

The ebola test data directory is also a Nextclade dataset (it contains a
pathogen.json declaring reference.fasta and genome_annotation.gff3). Sourcing
the annotation from it gives the same result as providing it explicitly with
--use-nextclade-gff-style.

  $ ${AUGUR} ancestral \
  >  --tree "$EBOLA/tree.nwk" \
  >  --alignment "$EBOLA/masked.fasta" \
  >  --nextclade-dataset "$EBOLA" \
  >  --genes GP L NP sGP ssGP VP24 VP30 VP35 VP40 \
  >  --translations "$EBOLA/translations/%GENE.fasta" \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  &> /dev/null

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$EBOLA/nt_muts.json" \
  >   ancestral_mutations.json
  {}

The path to the pathogen.json itself may be given instead of the directory.

  $ ${AUGUR} ancestral \
  >  --tree "$EBOLA/tree.nwk" \
  >  --alignment "$EBOLA/masked.fasta" \
  >  --nextclade-dataset "$EBOLA/pathogen.json" \
  >  --genes GP L NP sGP ssGP VP24 VP30 VP35 VP40 \
  >  --translations "$EBOLA/translations/%GENE.fasta" \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations_pathogen_json.json \
  >  &> /dev/null

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   ancestral_mutations.json \
  >   ancestral_mutations_pathogen_json.json
  {}

--nextclade-dataset and --annotation are mutually exclusive.

  $ ${AUGUR} ancestral \
  >  --tree "$EBOLA/tree.nwk" \
  >  --alignment "$EBOLA/masked.fasta" \
  >  --nextclade-dataset "$EBOLA" \
  >  --annotation "$EBOLA/genome_annotation.gff3" \
  >  --genes GP \
  >  --translations "$EBOLA/translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations.json
  ERROR: --nextclade-dataset and --annotation can not be used together.
  [2]

A dataset without a genome annotation can't be used for amino acid reconstruction.

  $ mkdir no-annotation
  $ cp "$EBOLA/reference.fasta" no-annotation/
  $ echo '{"files": {"reference": "reference.fasta"}}' > no-annotation/pathogen.json

  $ ${AUGUR} ancestral \
  >  --tree "$EBOLA/tree.nwk" \
  >  --alignment "$EBOLA/masked.fasta" \
  >  --nextclade-dataset no-annotation \
  >  --genes GP \
  >  --translations "$EBOLA/translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations.json
  ERROR: Nextclade dataset 'no-annotation' does not declare a genome annotation (files.genomeAnnotation).
  [2]
