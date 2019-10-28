# Sharing your analysis

[nextstrain.org](https://nextstrain.org) has a feature that allows you to share your own analysis through the nextstrain.org website.
This works using github as a repository for your analysis files.
To share an analysis, you need to create a repository.
The name of the repository will be what is used to access the results.
Within this repository, there should a folder called auspice which contains the output json files of the augur pipeline.
Importantly, the name of these files has to start with the name of the repository.
If so, you should be able to access your analysis via the nextstrain community feature.

As an example, lets look at one of our nextstrain community analysis.
The following link shows you an analysis we made a few month ago of many influenza B sequences:

[nextstrain.org/community/neherlab/allflu/B_ha](https://nextstrain.org/community/neherlab/allflu/B_ha)

The analysis files are hosted on our github page in the reposity

[github.com/neherlab/allflu](https://github.com/neherlab/allflu)

In this repository, you'll find the folder `auspice` which contains the files
```
allflu_B_ha_meta.json
allflu_B_ha_tree.json
```
Note that all files start with "allflu" which matches the name of the repository.
In fact, there are multiple analysis in this folder. The corresponding files all start with "allflu" but they differ in the viral lineage they correspond to:
```
allflu_B_ha_meta.json
allflu_B_ha_tree.json
allflu_h1n1_ha_meta.json
allflu_h1n1_ha_tree.json
allflu_h1n1pdm_ha_meta.json
allflu_h1n1pdm_ha_tree.json
allflu_h3n2_ha_meta.json
allflu_h3n2_ha_tree.json
```
All these can be accessed as

[nextstrain.org/community/neherlab/allflu/h1n1_ha](https://nextstrain.org/community/neherlab/allflu/h1n1_ha)

[nextstrain.org/community/neherlab/allflu/h3n2_ha](https://nextstrain.org/community/neherlab/allflu/h3n2_ha)

[nextstrain.org/community/neherlab/allflu/h1n1pdm_ha](https://nextstrain.org/community/neherlab/allflu/h1n1pdm_ha)


