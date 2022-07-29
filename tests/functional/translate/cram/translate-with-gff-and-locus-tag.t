Setup

  $ pushd "$TESTDIR" > /dev/null
  $ source _setup.sh

Translate amino acids for genes using a GFF3 file where the gene names are stored in a qualifier named "locus_tag".

  $ ${AUGUR} translate \
  >   --tree translate/data/tb/tree.nwk \
  >   --genes translate/data/tb/genes.txt \
  >   --vcf-reference translate/data/tb/ref.fasta \
  >   --ancestral-sequences translate/data/tb/nt_muts.vcf \
  >   --reference-sequence translate/data/tb/Mtb_H37Rv_NCBI_Annot.gff \
  >   --output-node-data $TMP/aa_muts.json \
  >   --alignment-output $TMP/translations.vcf \
  >   --vcf-reference-output $TMP/translations_reference.fasta
  Gene length of rrs_Rvnr01 is not a multiple of 3. will pad with N
  Read in 187 specified genes to translate.
  Read in 187 features from reference sequence file
  162 genes had no mutations and so have been be excluded.
  amino acid mutations written to .* (re)
  $ python3 "../../scripts/diff_jsons.py" translate/data/tb/aa_muts.json $TMP/aa_muts.json
  {}
