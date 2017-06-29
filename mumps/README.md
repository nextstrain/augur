## VIPR
* paramyxos -> genome search -> data to return = genome -> type "mumps" in search
* select complete genome only (n=104 as of july 2017)
* select all ->  download
* download genome fasta with custom format. Fields: GB-SN-SEG-DATE-HOST-COUN-SUBTY-VIRSPEC

## Clean FASTA
copy vipr fasta file to `fauna/data/mumps_vipr_full.fasta`
```
cd fauna
python mumps_restrict_vipr_fasta.py
```
The cleaned fasta is now available at `data/mumps_vipr.fasta`

## Upload to fauna
`python vdb/zika_upload.py -db vdb -v mumps --source genbank --locus genome --fname mumps_vipr.fasta`


## Download from fauna
`python vdb/zika_download.py -db vdb -v mumps --fstem mumps --resolve_method choose_genbank`


## AUGUR
```
cd augur/mumps
python prepare.mumps.py
python process.mumps.py
```
