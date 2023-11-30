Setup

  $ export AUGUR="${AUGUR:-$TESTDIR/../../../../bin/augur}"
  $ export DATA="$TESTDIR/../data"
  $ export SCRIPTS="$TESTDIR/../../../../scripts"

Translate amino acids for genes using a GFF3 file where the gene names are stored in a qualifier named "locus_tag".

  $ ${AUGUR} translate \
  >   --tree "${DATA}/tb/tree.nwk" \
  >   --genes "${DATA}/tb/genes.txt" \
  >   --vcf-reference "${DATA}/tb/ref.fasta" \
  >   --ancestral-sequences "${DATA}/tb/nt_muts.vcf" \
  >   --reference-sequence "${DATA}/tb/Mtb_H37Rv_NCBI_Annot.gff" \
  >   --output-node-data aa_muts.json \
  >   --alignment-output translations.vcf \
  >   --vcf-reference-output translations_reference.fasta
  Gene length of 'rrs' is not a multiple of 3. will pad with N
  Read in 187 specified genes to translate.
  Read in 187 features from reference sequence file
  162 genes had no mutations and so have been be excluded.
  amino acid mutations written to .* (re)

  $ python3 "${SCRIPTS}/diff_jsons.py" "${DATA}/tb/aa_muts.json" aa_muts.json \
  >  --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]"
  {}
