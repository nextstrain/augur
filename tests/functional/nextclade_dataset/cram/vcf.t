Setup

  $ source "$TESTDIR"/_setup.sh
  $ export ANC_DATA="$TESTDIR/../../ancestral/data/simple-genome"
  $ export DATA="$TESTDIR/../../translate/data/simple-genome"

Create a Nextclade dataset from the simple-genome reference and annotation.

  $ mkdir dataset
  $ cp "$ANC_DATA/reference.fasta" "$DATA/reference.gff" dataset/
  $ cat > dataset/pathogen.json <<~~
  > {"files": {"reference": "reference.fasta", "genomeAnnotation": "reference.gff"}}
  > ~~

For VCF input the dataset's reference is used as --vcf-reference. The genes of
the annotation don't have CDSs and so are used directly, giving the same result
as the default GFF parsing.

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$DATA/snps-inferred.vcf" \
  >  --nextclade-dataset dataset \
  >  --output-node-data aa_muts.json
  Read in 3 features from reference sequence file
  Validating schema of 'aa_muts.json'...
  amino acid mutations written to aa_muts.json

  $ python3 "$SCRIPTS/diff_jsons.py" \
  >   "$DATA/aa_muts.json" \
  >   aa_muts.json \
  >   --exclude-regex-paths "root\['annotations'\]\['.+'\]\['seqid'\]" "root['meta']['updated']"
  {}

--nextclade-dataset and --vcf-reference are mutually exclusive.

  $ ${AUGUR} translate \
  >  --tree "$ANC_DATA/tree.nwk" \
  >  --ancestral-sequences "$DATA/snps-inferred.vcf" \
  >  --nextclade-dataset dataset \
  >  --vcf-reference "$ANC_DATA/reference.fasta" \
  >  --output-node-data aa_muts.json
  ERROR: --nextclade-dataset and --vcf-reference can not be used together.
  [2]

The same applies to augur tree.

  $ ${AUGUR} tree \
  >  --alignment "$ANC_DATA/snps.vcf" \
  >  --nextclade-dataset dataset \
  >  --method fasttree \
  >  --output tree_dataset.nwk > /dev/null 2>&1

  $ ${AUGUR} tree \
  >  --alignment "$ANC_DATA/snps.vcf" \
  >  --vcf-reference "$ANC_DATA/reference.fasta" \
  >  --method fasttree \
  >  --output tree_explicit.nwk > /dev/null 2>&1

  $ diff tree_dataset.nwk tree_explicit.nwk
