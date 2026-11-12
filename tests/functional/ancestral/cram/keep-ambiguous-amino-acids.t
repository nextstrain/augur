Setup

  $ source "$TESTDIR"/_setup.sh

Mask the first 10 amino acids of one tip's ENV translation.

  $ cp "$TESTDIR"/../data/aa_sequences_PRO.fasta aa_sequences_PRO.fasta
  $ python3 -c '
  > from Bio import SeqIO
  > records = list(SeqIO.parse("'"$TESTDIR"'/../data/aa_sequences_ENV.fasta", "fasta"))
  > records[0].seq = "X" * 10 + records[0].seq[10:]
  > SeqIO.write(records, "aa_sequences_ENV.fasta", "fasta")
  > print(records[0].id)' > masked_tip.txt

Reconstruct amino acid sequences, explicitly requesting that ambiguous states are NOT inferred.
This works together with --output-translations and --report-inconsistent-translation.

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR"/../data/tree.nwk \
  >  --alignment "$TESTDIR"/../data/aligned.fasta \
  >  --annotation "$TESTDIR"/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations "aa_sequences_%GENE.fasta" \
  >  --keep-ambiguous \
  >  --report-inconsistent-translation \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations.json \
  >  --output-translations "ancestral_aa_sequences_%GENE.fasta" &> /dev/null

The masked tip keeps its ambiguous amino acids in the output.

  $ grep -A1 "^>$(cat masked_tip.txt)$" ancestral_aa_sequences_ENV.fasta | tail -1 | cut -c 1-12
  XXXXXXXXXX[^X]{2} (re)

  $ grep -A 2 "aa_sequences" ancestral_mutations.json
        "aa_sequences": {
          "ENV": .* (re)
          "PRO": .* (re)

With --infer-ambiguous (the default) the ambiguous amino acids are inferred.

  $ ${AUGUR} ancestral \
  >  --tree "$TESTDIR"/../data/tree.nwk \
  >  --alignment "$TESTDIR"/../data/aligned.fasta \
  >  --annotation "$TESTDIR"/../data/zika_outgroup.gb \
  >  --genes ENV PRO \
  >  --translations "aa_sequences_%GENE.fasta" \
  >  --seed 314159 \
  >  --output-node-data ancestral_mutations_inferred.json \
  >  --output-translations "inferred_aa_sequences_%GENE.fasta" &> /dev/null

  $ grep -A1 "^>$(cat masked_tip.txt)$" inferred_aa_sequences_ENV.fasta | tail -1 | cut -c 1-12 | grep -c X
  0
  [1]
