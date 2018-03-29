# Scores

Node-specific scores enable analysis and visual display of pathogenic characteristics such as the number of epitope mutations, age of infected individuals, or the local branching index (LBI).
We categorize scores into the following groups based on the data source from which they are calculated.

## Sequence scores

We calculate the following sequence scores from amino acid sequences.
We annotate each node of a given tree with a key in its `attr` attribute.

| Score | Attribute key | Description |
| epitope mutations | `ep` | Number of epitope mutations between the sequences of the root and the current node |
| non-epitope mutations | `ne` | Number of non-epitope mutations between the root and current node |
| receptor binding site mutations | `rb` | Number of mutations at experimentally-defined receptor binding sites |
| glycosylation sites | `glyc` | Number of glycosylation sites identified in the current node's sequence |

### Sequence masks

Mutation scores depend on previously defined sets of biologically-related amino acid positions in a given gene (e.g., epitopes or receptor binding sites in HA).
We define each set by a binary string that represents each position in the gene.
Positions included in a given set have a value of 1.
For example, H3N2 epitope mutations are defined by [Wolf. et al. 2006](https://www.ncbi.nlm.nih.gov/pubmed/17067369).
H3N2 receptor binding sites are defined by [Koel et al. 2013](https://www.ncbi.nlm.nih.gov/pubmed/24264991).

Masks are stored as a tab-delimited file in the metadata for a pathogen build (e.g., `builds/flu/metadata/ha_masks.tsv`).
The first column of this file specifies the name of the mask.
The second column specifies the binary string of the mask.

## Metadata scores

We annotate the following scores from metadata that accompanies each sequence record.

| Score | Attribute key | Description |
| age | `age` | Age of infected individual in years for tips or average age of individuals associated with tips for internal nodes |
| sex | `num_gender` | A numerical coding of the sex of the infected individual (male: -1, female: 1, and unknown: 0) |

## Phylogenetic scores

We annotate the following scores from the phylogeny rather than from individual node characteristics.

| Score | Attribute key | Description |
| local branching index | `lb` | Local branching index of a node with respect to a specific neighborhood on the phylogeny as defined by [Neher et al. 2013](https://elifesciences.org/articles/03568) |
