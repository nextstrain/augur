Setup

  $ export AUGUR="$TESTDIR/../../../bin/augur"
  $ export TB_DRM_DATA="$TESTDIR/../../builds/tb_drm/data"

A DRM file with the following columns does not exit with errors.

  $ cat >DRMs-valid.txt <<~~
  > GENE	SITE	ALT	DISPLAY_NAME	FEATURE
  > gyrB	461	N		Fluoroquinolones
  > gyrB	499	D		Fluoroquinolones
  > gyrA	88	C		Fluoroquinolones
  > ~~

  $ ${AUGUR} sequence-traits \
  >  --ancestral-sequences "$TB_DRM_DATA"/drm.vcf.gz \
  >  --vcf-reference "$TB_DRM_DATA"/ref.fasta \
  >  --features DRMs-valid.txt \
  >  --output-node-data output.json
  .* ParserWarning: Falling back to the 'python' engine because the 'c' engine does not support sep=None with delim_whitespace=False; you can avoid this warning by specifying engine='python'. (re)
  .* (re)
  This method may change in future! Please use 'augur sequence-traits -h' to check the latest options.
  Unfortunately this method currently only works with VCF input.
  sequence traits written to output.json

A DRM file with out the SITE column should produce a helpful error message.

  $ cat >DRMs-missing-site.txt <<~~
  > GENE	ALT	DISPLAY_NAME	FEATURE
  > gyrB	N		Fluoroquinolones
  > gyrB	D		Fluoroquinolones
  > gyrA	C		Fluoroquinolones
  > ~~

  $ ${AUGUR} sequence-traits \
  >  --ancestral-sequences "$TB_DRM_DATA"/drm.vcf.gz \
  >  --vcf-reference "$TB_DRM_DATA"/ref.fasta \
  >  --features DRMs-missing-site.txt \
  >  --output-node-data output.json
  .* ParserWarning: Falling back to the 'python' engine because the 'c' engine does not support sep=None with delim_whitespace=False; you can avoid this warning by specifying engine='python'. (re)
  .* (re)
  ERROR: DRMs-missing-site.txt has columns ['ALT', 'DISPLAY_NAME', 'FEATURE', 'GENE'], but expected ['ALT', 'FEATURE', 'GENE', 'SITE'] with optional ['DISPLAY_NAME'].
  This method may change in future! Please use 'augur sequence-traits -h' to check the latest options.
  Unfortunately this method currently only works with VCF input.
  [2]
