Setup

  $ source "$TESTDIR"/_setup.sh
  $ export DATA="$TESTDIR/../data/ebola"

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
  Reconstructing all 9 CDSs with translations: NP, VP35, VP40, GP, ssGP, sGP, VP30, VP24, L

  $ python3 "$TESTDIR/../../../../scripts/diff_jsons.py" \
  >   --exclude-regex-paths "\['seqid'\]" -- \
  >   "$DATA/nt_muts.json" \
  >   ancestral_mutations.json
  {}

  $ ls ancestral_*.fasta | wc -l | tr -d ' '
  9

CDSs without a (non-empty) translations file are skipped when --genes is omitted.
Tips without a sequence in a translations file (common for Nextclade
translations) are treated as fully ambiguous with a warning.

  $ mkdir translations
  $ cp "$DATA"/translations/*.fasta translations/
  $ rm translations/VP24.fasta
  $ : > translations/VP30.fasta
  $ python3 -c '
  > from Bio import SeqIO
  > records = list(SeqIO.parse("translations/GP.fasta", "fasta"))
  > SeqIO.write([r for r in records if r.id not in ("PP_00001QJ", "PP_00001C9")], "translations/GP.fasta", "fasta")' > /dev/null

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --translations "translations/{cds}.fasta" \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations_partial.json 2>&1 | grep "WARNING\|Reconstructing all"
  WARNING: Skipping 2 of 9 CDSs in '.*/genome_annotation.gff3' without a \(non-empty\) translations file: VP30, VP24 (re)
  WARNING: 2 of 10 tips have no sequence in 'translations/GP.fasta' and are treated as fully ambiguous for 'GP': PP_00001C9, PP_00001QJ
  Reconstructing all 7 CDSs with translations: NP, VP35, VP40, GP, ssGP, sGP, L

  $ python3 -c 'import json; print(sorted(json.load(open("ancestral_mutations_partial.json"))["nodes"]["PP_00001QJ"]["aa_muts"]))'
  ['GP', 'L', 'NP', 'VP35', 'VP40', 'sGP', 'ssGP']

A gene listed explicitly must have a translations file. This is checked before
any reconstruction.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --genes GP VP24 \
  >  --translations "translations/%GENE.fasta" \
  >  --output-node-data ancestral_mutations_error.json
  ERROR: The translations file 'translations/VP24.fasta' for gene 'VP24' does not exist.
  [2]

A gene listed explicitly must be in the annotation.

  $ ${AUGUR} ancestral \
  >  --tree "$DATA/tree.nwk" \
  >  --alignment "$DATA/masked.fasta" \
  >  --annotation "$DATA/genome_annotation.gff3" \
  >  --use-nextclade-gff-style \
  >  --genes gp \
  >  --translations "$DATA/translations/GP.fasta" \
  >  --output-node-data ancestral_mutations_error.json 2>&1 | grep ERROR
  ERROR: Gene\(s\) gp not found in '.*/genome_annotation.gff3'. Available genes/CDSs: NP, VP35, VP40, GP, ssGP, sGP, VP30, VP24, L (re)
  [2]
