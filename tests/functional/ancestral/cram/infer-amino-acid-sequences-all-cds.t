Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data/nextclade-gff"

If --genes is omitted, all CDSs in the annotation are reconstructed. The
translations may use Nextclade's {cds} placeholder instead of %GENE. This gives
the same result as listing all genes explicitly.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --translations "$DATA/translations/{cds}.fasta" \
  >  --infer-ambiguous \
  >  --inference joint \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  --output-translations "ancestral_{cds}.fasta" 2>&1 | grep "Reconstructing all"
  Reconstructing all 3 CDSs with translations: A, B, C

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$DATA/nt_muts.json" \
  >   ancestral_mutations.json
  {}

  $ ls ancestral_*.fasta
  ancestral_A.fasta
  ancestral_B.fasta
  ancestral_C.fasta

CDSs without a (non-empty) translations file are skipped when --genes is omitted.
Tips without a sequence in a translations file (common for Nextclade
translations) are treated as fully ambiguous with a warning.

  $ mkdir translations
  $ cp "$DATA"/translations/*.fasta translations/
  $ rm translations/C.fasta
  $ : > translations/B.fasta
  $ python3 -c '
  > from Bio import SeqIO
  > records = list(SeqIO.parse("translations/A.fasta", "fasta"))
  > SeqIO.write([r for r in records if r.id not in ("tip_1", "tip_6")], "translations/A.fasta", "fasta")' > /dev/null

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --translations "translations/{cds}.fasta" \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations_partial.json 2>&1 | grep "WARNING\|Reconstructing all"
  WARNING: Skipping 2 of 3 CDSs in '.*/genome_annotation.gff3' without a \(non-empty\) translations file: B, C (re)
  WARNING: 2 of 8 tips have no sequence in 'translations/A.fasta' and are treated as fully ambiguous for 'A': tip_1, tip_6
  Reconstructing all 1 CDSs with translations: A

  $ python3 -c 'import json; print(sorted(json.load(open("ancestral_mutations_partial.json"))["nodes"]["node_1"]["aa_muts"]))'
  ['A']

A gene listed explicitly must have a translations file. This is checked before
any reconstruction.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --genes A C \
  >  --translations "translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations_error.json
  ERROR: The translations file 'translations/C.fasta' for gene 'C' does not exist.
  [2]

A gene listed explicitly must be in the annotation.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --genes a \
  >  --translations "$DATA/translations/A.fasta" \
  >  --output-node-data ancestral_mutations_error.json 2>&1 | grep ERROR
  ERROR: Gene\(s\) a not found in '.*/genome_annotation.gff3'. Available genes/CDSs: A, B, C (re)
  [2]
