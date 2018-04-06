## Zika build

### How to run

* `cd` into this directory
* Make sure there's an (updated) fasta file in `./data`. Currently this is named `./data/full_dataset.singleline.aligned.pipeChar.fasta` and set in `wnv.prepare.py`.
* `python zika.prepare.py` (very quick) creates the JSON `prepared/WNV_NA.json`
* `python zika.process.py` (slow) creates the auspice-ready JSONs in `auspice`
* `cp auspice/WNV_NA_meta.JSON auspice/WNV_NA_tree.JSON ../../../auspice/data/` to test in auspice
* (to test in auspice) `cd ../../../auspice; npm run start:local`
* push to S3 via:
  * `cp auspice/WNV_NA_meta.JSON auspice/WNV_NA_tree.JSON ../../s3/`
  * `cd ../../s3`
  * `python ../scripts/s3.py push nextstrain-staging WNV_NA*` or `python ../scripts/s3.py push nextstrain-data WNV_NA*`
