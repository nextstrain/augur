# Zika tutorial -- fasta input

The tool-chain [*augur*](https://github.com/nextstrain/augur) is the bioinformatics engine of nextstrain and produces the files that can be visualized in the webbrowser using [*auspice*](https://github.com/nextstrain/auspice).
Augur consists of a number of tools that allow the user to filter and align sequences, build trees, and integrate the phylogenetic analysis with meta data.
The different tools are meant to be composable and the output of one tool will serve as the input of other tools.
We will work off the tutorial for Zika virus on the [nextstrain web site](https://nextstrain.org/docs/getting-started/zika-tutorial) and the github repository [nextstrain/zika-tutorial](https://github.com/nextstrain/zika-tutorial).

## Setup

To run this tutorial you'll need to [install augur](../guides/install/augur_install.md) and [install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

## Augur commands

As an example, we'll look that the `filter` command in greater detail.
This command allows you to selected various subsets of your input data for different types of analysis.
A simple example use of this command would be
```bash
augur filter --sequences data/sequences.fasta --metadata data/metadata.tsv --min-date 2012 --output filtered.fasta
```
This command will select all sequences with collection date in 2012 or later.
The filter command has a large number of options that allow flexible filtering for many common situations.
One such use-case is the exclusion of sequences that are known to be outliers (e.g.~because of sequencing errors, cell-culture adaptation, ...).
These can be specified in a separate file:
```
BRA/2016/FC_DQ75D1
COL/FLR_00034/2015
...
```
To drop such strains, you can pass the name of this file to the augur filter command:
```bash
augur filter --sequences data/sequences.fasta \
             --metadata data/metadata.tsv \
             --min-date 2012 \
             --exclude config/dropped_strains.txt \
             --output filtered.fasta
```
(To improve legibility, we have wrapped the command across multiple lines.)
If you run this command (you should be able to copy-paste this into your terminal), you should see that one of the sequences in the data set was dropped since its name was in the `drooped_strain.txt` file.

Another common filtering operation is subsetting of data to a achieve a more even spatio-temporal distribution or cut-down data set size to more manageable numbers.
The filter command allows you to select a specific number of sequences from specific groups, for example one sequence per month from each country:
```bash
augur filter \
  --sequences data/sequences.fasta \
  --metadata data/metadata.tsv \
  --min-date 2012 \
  --exclude config/dropped_strains.txt \
  --group-by country year month \
  --sequences-per-group 1 \
  --output filtered.fasta
```
This subsampling and filtering will reduce the number of sequences in this tutorial data set from 34 to 24.

## Chaining of augur commands with snakemake.

The output of one augur command serves as input to the next and hence these commands need to be executed in order.
This is a common pattern in bioinformatics and multiple tools -- so called workflow managers -- exist to facilitate this.
Within nextstrain, we have relied on [Snakemake](https://snakemake.readthedocs.io/en/stable/).

Snakemake breaks a workflow into a set of rules that are specified in a file called `Snakefile`.
Each rule takes a number of input files, specifies a few parameters, and produces output files.
A simple rule would look like this:
```python
rule filter:
    input:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv"
    output:
        sequences = "results/filtered.fasta"
    params:
        min_date = 2012
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output {output.sequences} \
            --min-date {params.min_date}
        """
```
This rule would produce `results/filtered.fasta` from the input files `data/sequences.fasta` and `data/metadata.tsv` using the `augur filter` command.
Note that we explicitly specify what is an input and what is an output file.
To filter our data, we would now call snakemake as
```bash
snakemake --cores 1 results/filtered.fasta
```
and snakemake will run the same command as specified above.

So far, this is just a complicated reformulation of what we did above, but the benefit of workflow management becomes clear once we add more steps.
The next natural step in our phylogenetic pipeline is aligning the filtered sequences and we define a rule `align`.
```bash
rule align:
    input:
        sequences = rules.filter.output.sequences,
        reference = "config/zika_outgroup.gb"
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment}
        """
```
If you now want to generate the alignment, you can type
```bash
snakemake --cores 1 results/aligned.fasta
```
and snakemake will

 * determine that `results/aligned.fasta` is an output of rule `align`
 * check whether all required input files are in place and run the necessary rules if not
 * run rule `align` and check whether the file `results/aligned.fasta` appeared.

If you supply a reference sequence to `augur align`, augur will include that reference sequence in the alignment and strip all insertions relative to that reference.
This will guarantee a consistent reference coordinate system and genome annotation later on.


### Construct the Phylogeny

Infer a phylogenetic tree from the multiple sequence alignment.
```bash
rule tree:
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """
```

The resulting tree is stored in [Newick format](http://evolution.genetics.washington.edu/phylip/newicktree.html).
Branch lengths in this tree measure nucleotide divergence.

### Get a Time-Resolved Tree
Most phylogenies you see on nextstrain are time-resolved, that is the branch lengths of the tree correspond to calendar time rather than evolutionary distance.
The default method to infer time-resolved phylogenies in nextstrain is [TreeTime](https://github.com/neherlab/treetime) and this analysis is an done using the augur command `refine` (this steps "refines" the existing tree...).
The corresponding rule for in the snakefile would look as follows:
```bash
rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = "data/metadata.tsv"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --timetree \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """
```
This command requires the output of the tree and align rules as input and will produce the file `results/tree.nwk` and `results/branch_lengths.json`.
The latter contains the inferred clock model and the inferred dates of all nodes, including internal nodes.
Each internal node has been given a unique name (e.g.~`NODE_0000001211`) which we will use in the following steps to attach more information to the tree.
The `refine` command has many different options that allow to specify clock rates, filter sequences that don't follow the molecular clock, calculate confidence intervals, etc.
We'll get into these details later.

## Annotate the Phylogeny

### Reconstruct Ancestral Traits

TreeTime can also infer ancestral traits from an existing phylogenetic tree and metadata annotating each tip of the tree.
The following command infers the region and country of all internal nodes from the time tree and original strain metadata.
As with the `refine` command, the resulting JSON output is indexed by strain or internal node name.

Specifying `--confidence` means that the confidence intervals for the reconstructed trait values will be estimated.

```bash
rule traits:
    input:
        tree = rules.refine.output.tree,
        metadata = "data/metadata.tsv"
    output:
        node_data = "results/traits.json",
    params:
        columns = "region country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """
```

### Infer Ancestral Sequences

Using `ancestral`, we can reconstruct what the ancesters of our samples' sequences must have looked like, and record the mutations that occurred on each branch.

Next, infer the ancestral sequence of each internal node and identify any nucleotide mutations on the branches leading to any node in the tree.

```bash
rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data}
        """
```

### Identify Amino-Acid Mutations

Identify amino acid mutations from the nucleotide mutations and a reference sequence with gene coordinate annotations.
The resulting JSON file contains amino acid mutations indexed by strain or internal node name and by gene name.
To export a FASTA file with the complete amino acid translations for each gene from each node’s sequence, specify the `--alignment-output` parameter in the form of `results/aligned_aa_%GENE.fasta`.

The reference sequence needs to be in Genbank format (if working with Fasta input) or GFF format (if working with VCF input) - see the Prerequisites page for more information on how to find and download a reference sequence from Genbank.

Open the Zika tutorial reference sequence `config/zika_outgroup.gb` in a text editor.
For translation, the important parts are the `source` and `CDS` sections that are part of the `Features` part of the file:

```bash
[...]
FEATURES             Location/Qualifiers
     source          1..10769
                     /collection_date="25-Oct-2013"
                     /country="French Polynesia"
                     /db_xref="taxon:64320"
                     /host="Homo sapiens"
                     /isolation_source="serum"
                     /mol_type="genomic RNA"
                     /organism="Zika virus"
                     /strain="PF13/251013-18"
[...]
     CDS             91..456
                     /product="capsid protein"
                     /gene="CA"
     CDS             457..735
                     /product="propeptide"
                     /gene="PRO"
     CDS             736..960
                     /product="membrane protein"
                     /gene="MP"
     CDS             961..2472
                     /product="envelope protein"
                     /gene="ENV"
[...]
```

The `source` section tells `augur` how long the whole genome is.
The `CDS` sections describe each gene (using `/gene=`) that should be translated by giving its name and its start and end location in the genome.

In your downloaded file from Genbank, you might notice that it has only one or no `CDS` sections, and that instead the sections you would like to translate as genes are described as `mat_peptide`, like the below:

```
     mat_peptide     91..465
                     /product="capsid protein"
     mat_peptide     466..735
                     /product="propeptide"
     mat_peptide     736..960
                     /product="membrane protein"
     mat_peptide     961..2475
                     /product="envelope protein"
```

In this case, you will have to manually modify this file to work in `augur`.
Anything you wish to translate as a gene should be changed from `mat_peptide` to `CDS` as its section name.
You will also need to add `/gene=` and a gene name to the sections you'd like to translate.
If you have one really long `CDS` section that covers most of the genome, we'd recommend not including this to translate (don't provide a `/gene=`) if you can include smaller genes instead, as such a long polyprotein is not very informative.

If you are working with VCF-input, the GFF file contains the same information, but `augur` will only translate those items coded as `gene`.

```bash
rule translate:
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = "config/zika_outgroup.gb"
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """
```

## Export the Results

Finally, collect all node annotations and metadata and export it all in auspice’s JSON format.
This command pulls together the metadata file and all of the output files generated by the previous steps, and combines it into two output files that contain the annotated tree (`"auspice/zika_tree.json"`) and associated metadata (`"auspice/zika_meta.json"`).

There is one input file below that is not metadata or output from a previous step: the config file (`"config/auspice_config.json"`).
This file is very important - it species the page title, maintainer(s), what filters are present, what places should be mapped by, and what the user is able to color the tree by ("color-bys").

Open up `"config/auspice_config.json"` in a text editor.
The `title` and `maintiner` fields should be easy to spot.
You can modify them if you wish (and should, if you are using your own data!).

The `filters` section specifies what traits will be available for users to filter the data by.
For example, if "country" is a filter, users will be able to only show data from the countries they select.

The `color_options` section specifies what users will be allowed to color the data by. You will almost always have the first two sections: `"gt"` and `"num_date"`, unless you are using trees without sequences or trees that aren't time-resolved.
After that, you can include any trait in your data.
If you have additional traits in your metadata, you will need to add them here to have them show up in auspice! You can copy the entry for `"country"` and re-name it to match the appropriate column from your metadata file.

The `geo` section informs auspice what traits should be used to draw the samples onto the map. Obviously, they must be location-related, like "country" and "region".
You can use other locations if that's more relevant to your data. We won't cover this right now, but can explain in detail as you come to this with your own data.

Don't worry about the `defaults` section for the moment.


The resulting tree and metadata JSON files are the inputs to the auspice visualization tool.

```bash
rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = "data/metadata.tsv",
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = "config/auspice_config.json"
    output:
        auspice_json = "auspice/zika.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """
```



