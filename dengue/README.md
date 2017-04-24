# DENGUE Pipeline Notes

`augur/dengue/dengue.py` can run either a serotype-specific build or an all-serotype build.  

**By default, `fauna/vdb/dengue_download.py` downloads all serotypes together.**  
`fauna$ python vdb/dengue_download.py -db vdb -v dengue`  
A single serotype can be downloaded like:   
`fauna$ python vdb/dengue_download.py -db vdb -v dengue --select serotype:1`

To run corresponding augur builds:  
`augur$ python dengue/dengue.py -s all` (default)  
`augur$ python dengue/dengue.py -s denv1`  

**N.B.:** By default, serotype-specific builds use the respective reference genomes specified by [LANL](https://hfv.lanl.gov/content/sequence/HFV/GenomeMapper/GenomeMapper.html).
The all-serotype build uses the reference genome from serotype 4.
