Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../data/simple-genome"


  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$DATA/snps-inferred.vcf" \
  >  --reference-sequence "$DATA/reference.gff" \
  >  --output-node-data "aa_muts.json" \
  >  --alignment-output aa_muts.vcf \
  >  --vcf-reference "$ANC_DATA/reference.fasta" \
  >  --vcf-reference-output reference.fasta
  Read in 3 features from reference sequence file
  Validating schema of 'aa_muts.json'...
  amino acid mutations written to aa_muts.json

  $ cat reference.fasta
  >gene1
  MPCG*
  >gene2
  MVK* (no-eol)

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   aa_muts.json \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {}

------------------------------   MISSING TEST   ----------------------------------
We should diff the `aa_muts.vcf` produced above against that in "${DATA}/aa_muts.vcf"
however there are a few subtleties:
- The translate-produced VCF uses diploid genotyping, rather than haploid
- The meta-lines are different (easy to skip via `grep -v '##'` or similar)
See <https://github.com/nextstrain/augur/issues/1356#issuecomment-1853194601> for
some additional context and discussion about the suitability of VCFs here.
----------------------------------------------------------------------------------