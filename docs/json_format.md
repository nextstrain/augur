# Format of augur output

Many augur commands output additional annotations and inferences made for nodes of the tree.
These data are stored as json files which serve as input to subsequent steps and the export step that combines all these annotations into a file for visualization.
The files output by different steps follow a common format:
```json
{
  "alignment": "results/aligned.fasta",
  "input_tree": "results/tree_raw.nwk",
  "clock": {
    "intercept": -2.2139560216719905,
    "rate": 0.0011001967850769103,
    "rtt_Tmrca": 2012.3272960820566
  },
  "nodes": {
    "1_0087_PF": {
      "branch_length": 0.0001694424712162532,
      "clock_length": 0.0001694424712162532,
      "date": "2013-12-31",
      "mutation_length": 0.00027877777449324505,
      "numdate": 2013.9993155382122,
      "raw_date": "2013-12-XX"
    },
    "1_0181_PF": {
      "branch_length": 0.00009101896681161861,
      "clock_length": 0.00009101896681161861,
      "date": "2013-12-31",
      "mutation_length": 0.0001858351915513057,
      "numdate": 2013.9993155382122,
      "raw_date": "2013-12-XX"
    },
    "NODE_0000030": {
      "branch_length": 0.0006468897126106284,
      "clock_length": 0.0006468897126106284,
      "date": "2016-06-06",
      "mutation_length": 0.0010229321826899364,
      "numdate": 2016.4318911009404
    }
  }
}
```
The json structure has one mandatory top-level field called `nodes` that contains inferences made for each node, possibly including internal nodes.
For that purpose, it is essential that internal nodes are consistently named.
The latter are assigned unique names starting with `NODE_` in the `refine` step and it is crucial that the topology and rooting of the tree remains untouched thereafter.
In addition to the mandatory `nodes` field, the json can contain additional information on the analysis performed, for example input parameters or other inferences such as the clock model.
