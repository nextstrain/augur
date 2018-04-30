## West Nile Virus (North America)

Please ensure the contents of this file are kept up to date with the static website (http://static.nextstrain.org/docs/01-builds/WNV-build/)

This build is not yet public, and the data has not been passed through sacra so there are some ad-hoc steps here.


### Datafiles needed:
* `data/full_dataset.singleline.aligned.pipeChar.fasta` (this is `cat data/full_dataset.singleline.aligned.fasta | tr '_' '|'` )
* `data/2018.04.08_WNV_headers_v2.csv`
* `data/state_list_v1.csv`


### Adding host / author information to the fasta header
`python enhance_fasta_headers.py` (this produces the file `data/WNV.fasta`)

### Generating the colour maps for different traits
`python make_colours.py` (this produces the file `colors.tsv`)

### Bioinformatics pipeline
* `python wnv.prepare.py -v 0` (very quick) creates the JSON `prepared/WNV_NA.json`
* `python wnv.process.py` (slow) creates the auspice-ready JSONs in `auspice`
* `cp auspice/WNV_NA_meta.JSON auspice/WNV_NA_tree.JSON ../../../auspice/data/` for testing in auspice
Then test in auspice (run `npm run start:local` in the auspice directory)


### Deploy to S3
* `cp auspice/WNV_NA_meta.JSON auspice/WNV_NA_tree.JSON ../../s3/`
* `cd ../../s3`
* `python ../scripts/s3.py push nextstrain-staging WNV_NA*` or `python ../scripts/s3.py push nextstrain-data WNV_NA*`
