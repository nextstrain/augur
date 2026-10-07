Setup

  $ source "$TESTDIR"/_setup.sh

Translating with the genome annotation sourced from a Nextclade dataset gives the
same result as providing it explicitly with --use-nextclade-gff-style.

  $ ${AUGUR} translate \
  >   --tree "$EBOLA/tree.nwk" \
  >   --ancestral-sequences "$EBOLA/nt_muts.json" \
  >   --nextclade-dataset "$EBOLA" \
  >   --output-node-data aa_muts_dataset.json > /dev/null 2>&1

  $ ${AUGUR} translate \
  >   --tree "$EBOLA/tree.nwk" \
  >   --ancestral-sequences "$EBOLA/nt_muts.json" \
  >   --reference-sequence "$EBOLA/genome_annotation.gff3" \
  >   --use-nextclade-gff-style \
  >   --output-node-data aa_muts_explicit.json > /dev/null 2>&1

  $ python3 "$SCRIPTS/diff_jsons.py" aa_muts_explicit.json aa_muts_dataset.json
  {}

  $ python3 -c 'import json; print(sorted(json.load(open("aa_muts_dataset.json"))["annotations"]))'
  ['GP', 'L', 'NP', 'VP24', 'VP30', 'VP35', 'VP40', 'nuc', 'sGP', 'ssGP']

One of --reference-sequence or --nextclade-dataset is required, but not both.

  $ ${AUGUR} translate \
  >   --tree "$EBOLA/tree.nwk" \
  >   --ancestral-sequences "$EBOLA/nt_muts.json" \
  >   --output-node-data aa_muts.json 2>&1 | tail -1
  augur translate: error: one of the arguments --reference-sequence --nextclade-dataset is required
  [2]

  $ ${AUGUR} translate \
  >   --tree "$EBOLA/tree.nwk" \
  >   --ancestral-sequences "$EBOLA/nt_muts.json" \
  >   --reference-sequence "$EBOLA/genome_annotation.gff3" \
  >   --nextclade-dataset "$EBOLA" \
  >   --output-node-data aa_muts.json 2>&1 | tail -1
  augur translate: error: argument --nextclade-dataset: not allowed with argument --reference-sequence
  [2]
