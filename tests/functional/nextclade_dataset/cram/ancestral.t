Setup

  $ source "$TESTDIR"/_setup.sh

The synthetic nextclade-gff test data directory is also a Nextclade dataset (it
contains a pathogen.json declaring reference.fasta and genome_annotation.gff3).
Sourcing the annotation from it gives the same result as providing it explicitly
with --use-nextclade-gff-style.

  $ ${AUGUR} ancestral \
  >  --tree "$DATASET/tree.nwk" \
  >  --alignment "$DATASET/masked.fasta" \
  >  --nextclade-dataset "$DATASET" \
  >  --genes A B C \
  >  --translations "$DATASET/translations/%GENE.fasta" \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  &> /dev/null

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$DATASET/nt_muts.json" \
  >   ancestral_mutations.json
  {}

The path to the pathogen.json itself may be given instead of the directory.

  $ ${AUGUR} ancestral \
  >  --tree "$DATASET/tree.nwk" \
  >  --alignment "$DATASET/masked.fasta" \
  >  --nextclade-dataset "$DATASET/pathogen.json" \
  >  --genes A B C \
  >  --translations "$DATASET/translations/%GENE.fasta" \
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
  >  --tree "$DATASET/tree.nwk" \
  >  --alignment "$DATASET/masked.fasta" \
  >  --nextclade-dataset "$DATASET" \
  >  --annotation "$DATASET/genome_annotation.gff3" \
  >  --genes A \
  >  --translations "$DATASET/translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations.json
  ERROR: --nextclade-dataset and --annotation can not be used together.
  [2]

A dataset without a genome annotation can't be used for amino acid reconstruction.

  $ mkdir no-annotation
  $ cp "$DATASET/reference.fasta" no-annotation/
  $ echo '{"files": {"reference": "reference.fasta"}}' > no-annotation/pathogen.json

  $ ${AUGUR} ancestral \
  >  --tree "$DATASET/tree.nwk" \
  >  --alignment "$DATASET/masked.fasta" \
  >  --nextclade-dataset no-annotation \
  >  --genes A \
  >  --translations "$DATASET/translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations.json
  ERROR: Nextclade dataset 'no-annotation' does not declare a genome annotation (files.genomeAnnotation).
  [2]
