## Download data from VIPR
* paramyxos -> genome search -> data to return = genome -> type "mumps" in search
* select complete genome only (n=104 as of july 2017)
* select all ->  download
* download genome fasta with custom format. Fields: GB-SN-SEG-DATE-HOST-COUN-SUBTY-VIRSPEC
* copy vipr fasta file to `fauna/data/mumps_vipr_full.fasta`

## Clean FASTA

```
cd fauna
python vdb/mumps_preprocess_fasta.py --vipr
```
The cleaned fasta is now available at `data/mumps_vipr.fasta`

_Note: this will eventually be incorporated into the Fauna upload script_

## Upload to fauna
`python vdb/mumps_upload.py -db vdb -v mumps --source genbank --locus genome --fname mumps_vipr.fasta`

## Update fauna database
This is needed to populate certain attributes such as author & paper title.
`python vdb/mumps_update.py -db vdb -v mumps --update_citations`

## Download from fauna
`python vdb/mumps_download.py -db vdb -v mumps --fstem mumps --resolve_method choose_genbank`


## AUGUR
```
cd augur/mumps
python prepare.mumps.py
python process.mumps.py
```
