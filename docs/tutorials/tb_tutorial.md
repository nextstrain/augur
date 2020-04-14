# MTb tutorials using vcf data

This tutorial explains how to create a Nextstrain build for Tuberculosis sequences.
However, much of it will be applicable to any run where you are starting with [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) files rather than [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.

As in the Zika fasta-input [tutorial](zika_tutorial), we'll build up a Snakefile step-by-step for each step of the analysis.

## Build steps
Nextstrain builds typically require the following steps:

* Preparing data
* Constructing a phylogeny
* Annotating the phylogeny
* Exporting the results

However, the commands that make up each step can vary depending on your pathogen and your analysis.
Here, we'll follow these steps:

 1. Prepare pathogen sequences and metadata
    1. Filter the sequences (remove unwanted sequences and/or sample sequences)
    2. Mask the sequences (exclude regions of the sequence that are unreliable)

 2. Construct a phylogeny
    1. Construct an initial tree (to get topology)
    2. Convert this into a time-resolved tree (to get accurate branch lengths)

 3. Annotate the phylogeny
    1. Infer ancestral sequences
    2. Translate genes and identify amino-acid changes
    3. Reconstruct ancestral states (like location or host)
    4. Identify clades on the tree
    5. Identify drug resistance mutations

 4. Export the final results into auspice-readable format

## Download Data

The data in this tutorial is public and is a subset of the data from Lee et al.'s 2015 paper [Population genomics of *Mycobacterium tuberculosis* in the Inuit](http://www.pnas.org/content/112/44/13609).
As location was anonymized in the paper, location data provided here was randomly chosen from the region for illustrative purposes.

First, download the Tuberculosis (TB) build which includes example data and a pathogen build script.
Then enter the directory you just cloned.

```bash
git clone https://github.com/nextstrain/tb.git
cd tb
conda activate nextstrain
```

## Prepare the Sequences

A Nextstrain build with VCF file input starts with:

* A VCF file containing all the sequences you want to include (variable sites only)
* A FASTA file of the reference sequence to which your VCF was mapped
* A tab-delimited metadata file

There are other files you will need if you want to perform certain steps, like masking.
If you are also working with TB sequences, you may be able to use the files provided here on your own data.
Otherwise, you'll need to provide files specific to your pathogen.

Here, our input file is compressed with gzip - you can see it ends with `.vcf.gz`.
However, `augur` can take gzipped or un-gzipped VCF files.
It can also produce either gzipped or un-gzipped VCF files as output (detected from the file ending you provide).
Here, we'll usually keep our output VCF files gzipped, by giving our output files endings like `.vcf.gz`, but you can specify `.vcf` instead.

All the data you need to make the TB build are in the `data` and `config` folders.

### Filter the Sequences

Sometimes you may want to exclude certain sequences from analysis.
You may also wish to downsample your data based on certain criteria. `filter` lets you do this.

For this example, we'll just exclude sequences in the file `dropped_strains.txt`.

We'll need to specify these starting files at the top of our Snakefile:
```bash
seq_file = "data/lee_2015.vcf.gz",
meta_file = "data/meta.tsv",
exclude_file = "config/dropped_strains.txt"
```

And we'll add this as our first rule:
```bash
rule filter:
    input:
        seq = seq_file,
        meta = meta_file,
        exclude = exclude_file
    output:
        "results/filtered.vcf.gz"
    shell:
        """
        augur filter --sequences {input.seq} \
            --metadata {input.meta} \
            --exclude {input.exclude} \
            --output {output}
        """
```

Now run filter. If you are using the Snakefile included with the TB tutorial, you can run:
```bash
snakemake filter
```

If you have created your own Snakefile, you'll need to specify its name. For example, if it is called `TB_snakefile`, you would run:
```bash
snakemake -s TB_snakefile filter
```

### Mask the Sequences

There may be regions in your pathogen sequences that are unreliable.
For example, areas that are hard to map because of repeat regions.
Often, these are excluded from analysis so that incorrect calls in these areas don't influence the results.
The areas to be masked are specified in a BED-format file.
This is a standard, tab-delimited format with five columns: Chrom, ChomStart, ChromEnd, locus tag, and Comment.
You can open up `config/Locus_to_exclude_Mtb.bed` in the TB tutorial to see the file format.

The first, fourth, and fifth columns (Chrom, locus tag, and Comment) can be blank or contain anything - they will be ignored.
All sites between each ChromStart and ChromEnd will be removed from the analysis.

We'll need to add this BED-format file to the top of the Snakefile (below the files already there):
```bash
mask_file = "config/Locus_to_exclude_Mtb.bed"
```

Now we can add the `mask` rule:
```bash
rule mask:
    input:
        seq = rules.filter.output,
        mask = mask_file
    output:
       "results/masked.vcf.gz"
    shell:
        """
        augur mask --sequences {input.seq} \
            --mask {input.mask} \
            --output {output}
        """
```

## Construct the Phylogeny

Now our sequences are ready to start analysis.

With VCF files, we'll do this in two steps that are slightly different from FASTA-input.
1. First, we'll use only the variable sites to construct a tree quickly. This will give us the topology, but the branch lengths will be incorrect.
2. Next, we'll consider the entire sequence to correct our branch lengths.
At the same time, the sample date information will be used to create a time-resolved tree.

### Get the Topology

You can use different tree-building programs to build your initial tree, and specify some parameters.
Here, we'll use IQTree.
We specify it here with the argument `--method`, but it's also the default.

In `tree`, we pass in the VCF file and the reference it was mapped to.
We also pass in a list of sites that we'd like to exclude from building the topology (optional).
These are sites associated with drug-resistance mutations that can influence the topology.
We exclude them here, but they'll be allowed to influence branch length and be included in ancestral sequence reconstruction later.

We must add the reference sequence our VCF file was mapped to, and our list of sites to exclude from tree-building to the top of the Snakefile:
```bash
ref_file = "data/ref.fasta"
sites_file = "config/drm_sites.txt"
```

And add the `tree` rule to the Snakefile:
```bash
rule tree:
    input:
        aln = rules.mask.output,
        ref = ref_file,
        sites = sites_file
    output:
        "results/tree_raw.nwk"
    params:
        method = 'iqtree'
    shell:
        """
        augur tree --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --method {params.method} \
            --exclude-sites {input.sites} \
            --output {output}
        """
```

### Fix Branch Lengths & Get a Time-Resolved Tree

Now we'll use the topology from `tree`, but get more accurate branch lengths and a time-resolved tree.
This adjusts branch lengths in the tree to position tips by their sample date and infer the most likely time of their ancestors, using [TreeTime](https://github.com/neherlab/treetime).
There are many options that can be specified here in `refine` to help you get a good tree.

`refine` will produce as output:
* another tree (newick format)
* a JSON format file with the inferred dates and mutations on each node/branch

```bash
rule refine:
    input:
        tree = rules.tree.output,
        aln = rules.mask.output,
        metadata = meta_file,
        ref = ref_file
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json",
    params:
        root = 'min_dev',
        coal = 'opt'
    shell:
        """
        augur refine --tree {input.tree} \
            --alignment {input.aln} \
            --vcf-reference {input.ref} \
            --metadata {input.metadata} \
            --timetree \
            --root {params.root} \
            --coalescent {params.coal} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data}
        """
```

In addition to assigning times to internal nodes, the `refine` command filters tips that are likely outliers.
Branch lengths in the resulting Newick tree measure adjusted nucleotide divergence.
All other data inferred by TreeTime is stored by strain or internal node name in the JSON file.

## Annotate the Phylogeny

Now that we have an accurate tree and some information about the ancestral sequences, we can annotate some interesting data onto our phylogeny.
TreeTime can infer ancestral sequences and ancestral traits from an existing phylogenetic tree and metadata to annotate each tip of the tree.

### Infer Ancestral Sequences

We can reconstruct the ancestral sequences for the internal nodes on our phylogeny and identify any nucleotide mutations on the branches leading to any node in the tree.

For VCF runs, `ancestral` will produce another VCF that contains the reconstructed sequence of all the internal nodes and the sequences from the tip nodes, as well as a JSON-format file that contains nucleotide mutation information for each node.

```bash
rule ancestral:
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output,
        ref = ref_file
    output:
        nt_data = "results/nt_muts.json",
        vcf_out = "results/nt_muts.vcf"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral --tree {input.tree} \
            --alignment {input.alignment} \
            --vcf-reference {input.ref} \
            --inference {params.inference} \
            --output-node-data {output.nt_data} \
            --output-vcf {output.vcf_out}
        """
```

### Identify Amino-Acid Mutations

With `translate` we can identify amino acid mutations from the nucleotide mutations and a GFF file with gene coordinate annotations.
The resulting JSON file contains amino acid mutations indexed by strain or internal node name and by gene name.
`translate` will also produce a VCF-style file with the amino acid changes for each gene and each sequence, and FASTA file with the translated 'reference' genes which the VCF-style file 'maps' to.

Because of the number of genes in TB, we will only translate genes associated with drug resistance to save time.
We can pass in a list of genes to translate using `--genes`.
Note that the `--reference-sequence` option is how you pass in the GFF file with the gene coordinates.

We'll need to add the GFF file with the gene annotations and the file with a list of genes to translate to the list of files at the top of the Snakefile:
```bash
generef_file = "config/Mtb_H37Rv_NCBI_Annot.gff",
genes_file = "config/genes.txt"
```

```bash
rule translate:
    input:
        tree = rules.refine.output.tree,
        ref = ref_file,
        gene_ref = generef_file,
        vcf = rules.ancestral.output.vcf_out,
        genes = genes_file
    output:
        aa_data = "results/aa_muts.json",
        vcf_out = "results/translations.vcf",
        vcf_ref = "results/translations_reference.fasta"
    shell:
        """
        augur translate --tree {input.tree} \
            --vcf-reference {input.ref} \
            --ancestral-sequences {input.vcf} \
            --genes {input.genes} \
            --reference-sequence {input.gene_ref} \
            --output-node-data {output.aa_data} \
            --alignment-output {output.vcf_out} \
            --vcf-reference-output {output.vcf_ref}
        """
```

### Reconstruct Ancestral States

`traits` can reconstruct the probable ancestral state of traits like location and host (or others).
This is done by specifying a column or columns in the meta-data file.

`--confidence` will give confidence estimates for the reconstructed states.
The output will be a JSON file with the state (and confidence, if specified) information for each node.

```bash
rule traits:
    input:
        tree = rules.refine.output.tree,
        meta = meta_file
    output:
        "results/traits.json"
    params:
        traits = 'location'
    shell:
        """
        augur traits --tree {input.tree} \
            --metadata {input.meta} \
            --columns {params.traits} \
            --output-node-data {output}
        """
```

### Identify Specified Clades

In the [original paper](http://www.pnas.org/content/112/44/13609), the authors identified 'sublineages' within the dataset.
We can add these to our dataset as 'clades' by defining the sublineages with amino-acid or nucleotide mutations specific to that sublineage, given here in the file `config/clades.tsv`.
Open it up in a text editor to have a look at the format.

The `clades.tsv` file must be tab-delimited with four columns: clade, gene, site, and alt.
The 'clade' column gives the name of the clade being defined - you can have more than one row per clade - it will only be defined from the branch where all criteria are met.
The 'gene' and 'site' columns specify the gene (or `nuc` for nucleotide) and location (by AA position in the gene, or nucleotide position in the genome) where the branch must have the 'alt' (4th column) value to be considered this clade.

As clades, these sublineages will be labelled and we'll be able to color the tree by them.

You can specify clades for your own data by first doing a run without clades, then mousing over branches where you'd like to start defining a clade to see what mutations are present.

We'll need to add the file that defines the clades to the top of our Snakefile:
```bash
clades_file = "config/clades.tsv"
```

```bash
rule clades:
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.aa_data,
        nuc_muts = rules.ancestral.output.nt_data,
        clades = clades_file
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

### Identify Drug Resistance Mutations

`sequence-traits` can identify any trait associated with particular nucleotide or amino-acid mutations, not just drug resistance mutations (DRMs).

This dataset doesn't actually contain any drug resistance mutations, but identifying such mutations is often of interest to those working on tuberculosis.
Here, we'll run this step as an example, even though it won't add anything to the tree for this dataset.

Open up the `config/DRMs-AAnuc.tsv` file to see the format of a file that specifies sequence traits.
It contains five columns: GENE, SITE, ALT, DISPLAY_NAME, and FEATURE. DISPLAY_NAME can be blank.

For drug resistance, we list the gene, the AA position in the gene, the AA mutation that confers resistance (you can list a site multiple times if multiple bases give resistance), and the name of the drug this mutation gives resistance to:
```bash
GENE	SITE	ALT	DISPLAY_NAME	FEATURE
gyrB	461	N		Fluoroquinolones
gyrB	499	D		Fluoroquinolones
rpoB	432	E		Rifampicin
rpoB	432	K		Rifampicin
```
We can leave DISPLAY_NAME blank, as auspice will by default display the gene, site, and original and alternative base.

For mutations outside of protein-coding genes, we can specify their position using nucleotides:
```bash
GENE	SITE	ALT	DISPLAY_NAME	FEATURE
nuc	1472749	A	rrs: C904A	Streptomycin
nuc	1473246	G	rrs: A1401G	Amikacin Capreomycin Kanamycin
nuc	1673423	T	fabG1: G-17T	Isoniazid Ethionamide
nuc	1673425	T	fabG1: C-15T	Isoniazid Ethionamide
```
In the literature, these mutations are still referred to by their position within non-protein-coding genes (`rrs`) or location near genes (`-17 fabG1`), not their nucleotide location.
We can ensure auspice displays the more useful common nomenclature by giving entries for the DISPLAY_NAME column.

`sequence-traits` will return a value for each "feature" - for example, all the mutations on the tree that lead to resistance to Streptomycin.
It will also generate a count either of the total number of "features" each node has (ex: the total number of drugs a sequence is resistant to), or the total number or mutations specified in the file each node has (ex: the total number of DRMs a sequence has, even if some are for the same drug).
You can specify a name for this count using the `--label` argument (here: "Drug_Resistance").
The `--count` argument value specifies whether to count the number of traits (ex: drugs resistant to) (use `traits`) or number of overall mutations (use `mutations`).

We'll need to add the file that defines the sequence traits (DRMs) to the top of our Snakefile:
```bash
drms_file = "config/DRMs-AAnuc.tsv"
```

```bash
rule seqtraits:
    input:
        align = rules.ancestral.output.vcf_out,
        ref = ref_file,
        trans_align = rules.translate.output.vcf_out,
        trans_ref = rules.translate.output.vcf_ref,
        drms = drms_file
    output:
        drm_data = "results/drms.json"
    params:
        count = "traits",
        label = "Drug_Resistance"
    shell:
        """
        augur sequence-traits \
            --ancestral-sequences {input.align} \
            --vcf-reference {input.ref} \
            --translations {input.trans_align} \
            --vcf-translate-reference {input.trans_ref} \
            --features {input.drms} \
            --count {params.count} \
            --label {params.label} \
            --output-node-data {output.drm_data}
        """
```

## Export the Results

Finally, collect all node annotations and metadata and export it all in auspice’s JSON format.
The resulting tree and metadata JSON files are the inputs to the auspice visualization tool.

The names of the output tree and meta data files are here specified by a rule called `all` at the beginning of our Snakefile.
It should be even before the list of files, and looks like this:
```bash
rule all:
    input:
        auspice_tree = "auspice/tb_tree.json",
        auspice_meta = "auspice/tb_meta.json"
```

This rule tells Snakemake what the final output of our entire run should look like.
It will run all rules necessary to produce these files, so they should be the names of your final step.
If you have an "all" rule, you can run your entire analysis just by running `snakemake` or `snakemake --snakefile Snakefile2` (if the name of your Snakefile is not 'Snakefile').

We'll need to add a few remaining files to our list of files at the start of our Snakefile:
```bash
colors_file = "config/color.tsv",
config_file = "config/config.json",
geo_info_file = "config/lat_longs.tsv"
```

The `color.tsv` file is optional, but allows us to specify our own colors for particular traits.
If you open it up, you can see that we choose our own colors for values in 'region', 'country', 'location' and 'clade_membership'.
If you don't supply a `color.tsv` file, auspice will choose colors for you.
This can be the simplest way to start - then you can add colors for any traits where you don't like what auspice has chosen.

The `lat_longs.tsv` file contains the latitudes and longitudes for the geographic locations of your data, and may or may not be needed for your data.
Augur contains many latitudes and longitudes for countries and regions, but if you want to specify data at a different level (state, province, county, city), you can include your own file as well (it will be used in addition to the defaults, so country location can still be retrieved from the augur file, for example).
At the bottom of the `config/lat_longs.tsv` file in the TB tutorial, notice there are entries for 'location', listing each village.

The `config.json` file should be familiar from the Zika Fasta-input [tutorial](zika_tutorial).
There are a couple of changes worth pointing out, though.

Since all the samples come from the region of North America and the country of Canada, we don't include these anywhere in our data - all the samples would be the same.
Instead, we have 'location' as a `color_options` entry, and also as our `geo` (where the samples will be drawn on the map), and as a `filters` option.

We also have a `color_options` entry for 'clade_membership', since we designated clades with the `clades` rule.
The trait is added to our tree as `clade_membership` which is why this is the name of the option and the `key` value, but we could set the `legendTitle` and `menuItem` to be anything we wish, if we wanted.


```bash
rule export:
    input:
        tree = rules.refine.output.tree,
        metadata = meta_file,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output,
        nt_muts = rules.ancestral.output.nt_data,
        aa_muts = rules.translate.output.aa_data,
        drms = rules.seqtraits.output.drm_data,
        color_defs = "config/colors.tsv",
        config = "config/config.json",
        geo_info = "config/lat_longs.tsv",
        clades = rules.clades.output.clade_data
    output:
        auspice_json = "auspice/tb.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.drms} {input.aa_muts} {input.nt_muts} {input.clades} \
            --auspice-config {input.config} \
            --colors {input.color_defs} \
            --lat-longs {input.geo_info} \
            --output {output.auspice_json} \
            """
```


As mentioned previously, this dataset has no drug resistance, so it's not included in the `config.json` file to display, even though we ran the `sequence-traits` rule.
If you did have drug resistance information that you wanted to display, you would need to add it to the `config.json` file as `color_options`.

First, you would want to add a color-by for the total number of drugs each node is resistant to.
Since we gave the label 'Drug_Resistance' when we ran the rule, this will be the name of the option, and the `key`, but we can make the `menuItem` and `legendTitle` different if we wish:
```
  "Drug_Resistance": {
   "menuItem": "Drug_Resistance",
   "legendTitle": "Drug Resistance",
   "type": "discrete",
   "key": "Drug_Resistance"
  },
```
If you had given a different label when you ran the rule, you would change this entry to match.

You would then need an option for each drug where you have resistance information (or each FEATURE where you have information).
For example, to show the mutations present that confer resistance to Streptomycin and Rifampicin:
```
  "Streptomycin": {
   "menuItem": "Streptomycin",
   "legendTitle": "Streptomycin Resistance",
   "type": "discrete",
   "key": "Streptomycin"
  },
  "Rifampicin": {
   "menuItem": "Rifampicin",
   "legendTitle": "Rifampicin Resistance",
   "type": "discrete",
   "key": "Rifampicin"
  },
```
You would need an entry for every FEATURE in your original file (though you could then remove any that had no information on the tree).

