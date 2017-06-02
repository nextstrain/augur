# seasonal flu (h3n2 // h1n1pdm // vic // yam)

### status
* basic prepare + process working for h3n2 / ha only.
* single pass (no repeated) subsampling including titer information
* no titers
* no frequencies


### summary of the necessary files

| File         | Details           |
| ------------- | ------------- |
| `flu.prepare.py`    | Prepare script.      |
| `flu.process.py`    | Process script.      |
| `flu_info.py`      | Holds (a lot of) information about sequences to drop, reference genomes etc. Used by prepare.  |
| `../../fauna/flu_<SERO>_<SEG>.fasta` | (Fauna) fasta file      |
| `colors.flu.tsv` | color maps      |


### how to run
* download fauna files something like this:
```
cd fauna
SERO="h3n2"
SEG="ha"
python vdb/flu_download.py -db vdb -v flu --select locus:${SEG} lineage:seasonal_${SERO} --fstem ${SERO}_${SEG}
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


### To Do:

* process script needs these calls after the geo inference:
```
flu.estimate_tree_frequencies(pivots=pivots)
for region in regions:
    flu.estimate_tree_frequencies(region=region)
flu.dump()

flu.HI_model(criterium = lambda x:len(x.aa_mutations['HA1'])>0)
H3N2_scores(flu.tree.tree)
flu.dump()
flu.matchClades(clade_designations[params.lineage])
flu.export(extra_attr=['serum'], ...)
flu.HI_export()

# fn defined in seasonal_flu.py
plot_sequence_count(flu, store_data_path+'tip_count.png')
```
