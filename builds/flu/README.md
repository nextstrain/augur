# seasonal flu (h3n2 // h1n1pdm // vic // yam)

### status
* all lineages working
* references only defined for HA (so far)
* subsampling is always a single pass and cannot be repeated
* titer models are broken
* H3N2_scores is not implemented
* tree frequencies have not been tested properly
* matchClades has not been properly tested


### summary of the necessary files

| File         | Details           |
| ------------- | ------------- |
| `flu.prepare.py`    | Prepare script.      |
| `flu.process.py`    | Process script.      |
| `flu_info.py`      | Holds (a lot of) information about sequences to drop, reference genomes etc. Used by prepare.  |
| `../../fauna/flu_<SERO>_<SEG>.fasta` | (Fauna) fasta file      |
| `colors.flu.tsv` | color maps      |


### how to run
* download fauna files like this with subtype `h3n2`, `h1n1pdm`, `vic` or `yam` and segment `ha`, `na`, etc...:
```
cd fauna
python vdb/flu_download.py -db vdb -v flu --select lineage:seasonal_h3n2 locus:ha --fstem h3n2_ha
```

* prepare fauna fasta -> JSON ready to be analysed
```
python flu.prepare.py [options]
```

* process (prepared) JSON -> intermediate files + auspice JSONs
```
python flu.process.py [options]
cp auspice/* ../../auspice/data/
```
