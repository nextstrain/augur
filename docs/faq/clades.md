# Labeling `clades`

Clades in phylogenetic trees are often named to facilitate discussion of genetic diversity, see for example [seasonal influenza on nextstrain](https://nextstrain.org/flu).
Augur has a command to determine the position of such clade labels and assign sequences to clades.
The definition of these clades are provided in a tab-delimited file (tsv) using the following format:
```
clade	gene	site	alt
3b	HA1	145	S
3b	HA1	312	S
3b	nuc	1671	G
3c	HA1	48	I
3c	HA1	45	N
3c	nuc	456	T
3c2	HA2	160	N
3c2	nuc	693	A
```
Each line specifies a sequence feature of a clade.
The first column specified the name of the clade to which the feature belongs, the column `gene` can be any annotated gene or the underlying nucleotide sequence (`nuc`), the column `site` specifies the position (numbering starting at 1), while the column `alt` specified the state.
A clade if often defined by multiple such sequence features.

The augur command `clades` can be used to annotate such clades in your tree and a rule in a Snakefile would look this
```bash
rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        clades = "config/clades.tsv"
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output-node-data {output.clade_data}
        """
```
As input, this command requires the tree, the output of the ancestral reconstruction steps and the translation step (assuming your clades are defined using translations), as well as the file with clade definitions.

The output of this command is a json file in the common augur format that specifies `clade_membership` for each node in the tree.
Nodes that didn't match any clade definition will be left `unassigned`.
Internal nodes that form the root of clades will have an additional field `clade_annotation` that auspice uses to label branches in the tree.


