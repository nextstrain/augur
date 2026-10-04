Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data/nextclade-gff"

Infer ancestral nucleotide and amino acid sequences using a Nextclade-style GFF
of a synthetic 99 nt genome. CDS A is nested under a gene and an mRNA, CDS B is
split over two rows sharing an ID (with a codon spanning the junction) which
are joined into a single feature, and CDS C is on the minus strand.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --genes A B C \
  >  --use-nextclade-gff-style \
  >  --translations "$DATA/translations/%GENE.fasta" \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  &> /dev/null

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$DATA/nt_muts.json" \
  >   ancestral_mutations.json
  {}

  $ python3 -c '
  > import json
  > d = json.load(open("ancestral_mutations.json"))
  > for gene in "ABC":
  >     print(gene, {k: v for k, v in d["annotations"][gene].items() if k != "seqid"})
  > for name, node in sorted(d["nodes"].items()):
  >     print(name, node["muts"], {gene: muts for gene, muts in node["aa_muts"].items() if muts})'
  A {'end': 27, 'start': 4, 'strand': '+', 'type': 'CDS'}
  B {'segments': [{'end': 44, 'start': 31}, {'end': 63, 'start': 51}], 'strand': '+', 'type': 'CDS'}
  C {'end': 93, 'start': 67, 'strand': '-', 'type': 'CDS'}
  node_1 ['G10C', 'A99G'] {'A': ['E3Q']}
  node_2 ['G47C', 'A51G'] {'B': ['I5M']}
  node_3 ['G80C'] {'C': ['T5S']}
  node_4 ['A2G', 'A20G'] {'A': ['Q6R']}
  root [] {}
  tip_1 ['A35G'] {'B': ['Q2R']}
  tip_2 [] {}
  tip_3 ['C46-', 'C47-', 'A48-'] {}
  tip_4 [] {}
  tip_5 ['G70C'] {'C': ['D8E']}
  tip_6 ['T57C'] {}
  tip_7 [] {}
  tip_8 [] {}

The Nextclade GFF parsing flag requires an annotation file.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --use-nextclade-gff-style \
  >  --output-node-data ancestral_mutations.json
  ERROR: --use-nextclade-gff-style requires an --annotation file.
  [2]
