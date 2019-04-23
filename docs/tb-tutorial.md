# TB Tutorial

This tutorial explains how to create a Nextstrain build for Tuberculosis sequences.
However, much of it will be applicable to any run where you are starting with [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) files rather than [FASTA](https://en.wikipedia.org/wiki/FASTA_format) files.

We will first make the build step-by-step using an example data set.
Then we will see how to automate this stepwise process by defining a pathogen build script which contains the commands we will run below.

If you haven't already worked through the [Quickstart](quickstart), you may want to back up and begin there before continuing with this tutorial.

Note that we will not use all the bioinformatics commands possible with Nextstrain.
After running this tutorial, you may want to read more about the [bioinformatics commands offered by Nextstrain](/docs/bioinformatics).

## Setup

If you have the [Nextstrain command-line interface (CLI) tool](https://github.com/nextstrain/cli) installed from following the [Quickstart](quickstart), then you're all set!

Otherwise, you'll need to either [install the CLI](quickstart#set-up-your-computer) or [install the Nextstrain components](../getting-started/installation) individually.

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

```
git clone https://github.com/nextstrain/tb.git
cd tb
```

Next, if you're using the Nextstrain CLI tool, use it to enter the Nextstrain build environment by running:

```
nextstrain shell .
```

Note the dot (`.`) as the last argument; it is important and indicates that your current directory (`tb/`) is the build directory.
Your command prompt will change to indicate you are in the build environment.
(If you want to leave the build environment, run the command `exit`.)

If instead you installed Nextstrain components [using Conda](../getting-started/installation#install-augur-with-conda), remember to activate your environment by running:

```
conda activate nextstrain
```

## Prepare the Sequences

A Nextstrain build with VCF file input starts with:
* A VCF file containing all the sequences you want to include (variable sites only)
* A FASTA file of the reference sequence to which your VCF was mapped
* A tab-delimited metadata file _we need better info about what format this should be..._

There are other files you will need if you want to perform certain steps, like masking.
If you are also working with TB sequences, you may be able to use the files provided here on your own data.
Otherwise, you'll need to provide files specific to your pathogen.

Here, our input file is compressed with gzip - you can see it ends with `.vcf.gz`.
However, `augur` can take gzipped or un-gzipped VCF files.
It can also produce either gzipped or un-gzipped VCF files as output.
Here, we'll usually keep our VCF files gzipped, by giving our output files endings like `.vcf.gz`, but you can specify `.vcf` instead.

All the data you need to make the TB build is in the `data` and `config` folders.

### Filter the Sequences

Sometimes you may want to exclude certain sequences from analysis.
You may also wish to downsample your data based on certain criteria. `filter` lets you do this.
For more information on all the ways you can filter data, see [the bioinformatics commands reference](/docs/bioinformatics/commands#filter).

For this example, we'll just exclude sequences in the file `dropped_strains.txt`.

First, we'll make a folder to hold all of our results from each step:

```
mkdir -p results
```

Now run filter:

```
augur filter \
    --sequences data/lee_2015.vcf.gz \
    --metadata data/meta.tsv \
    --exclude config/dropped_strains.txt \
    --output results/filtered.vcf.gz
```

### Mask the Sequences

There may be regions in your pathogen sequences that are unreliable.
For example, areas that are hard to map because of repeat regions.
Often, these are excluded from analysis so that incorrect calls in these areas don't influence the results.
The areas to be masked are specified in a BED-format file.

```
augur mask \
    --sequences results/filtered.vcf.gz \
    --mask config/Locus_to_exclude_Mtb.bed \
    --output results/masked.vcf.gz
```

## Construct the Phylogeny

Now our sequences are ready to start analysis.

With VCF files, we'll do this in two steps that are slightly different from FASTA-input.
1. First, we'll use only the variable sites to construct a tree quickly. This will give us the topology, but the branch lengths will be incorrect.
2. Next, we'll consider the entire sequence to correct our branch lengths.
At the same time, the sample date information will be used to create a time-resolved tree.

### Get the Topology

You can use different tree-building programs to build your initial tree, and specify some parameters.
You can see all the options for `tree` in [the bioinformatics commands reference](/docs/bioinformatics/commands#tree).

Here, we pass in the VCF file and the reference it was mapped to.
We also pass in a list of sites that we'd like to exclude from building the topology.
These are sites associated with drug-resistance mutations that can influence the topology.
We exclude them here, but they'll be allowed to influence branch length and be included in ancestral sequence reconstruction later.
Finally, we use `iqtree` as the method to build the tree here.

```
augur tree \
    --alignment results/masked.vcf.gz \
    --vcf-reference data/ref.fasta \
    --exclude-sites config/drm_sites.txt \
    --method iqtree \
    --output results/tree_raw.nwk
```

### Fix Branch Lengths & Get a Time-Resolved Tree

Now we'll use the topology from `tree`, but get more accurate branch lengths and a time-resolved tree.
This adjusts branch lengths in the tree to position tips by their sample date and infer the most likely time of their ancestors, using [TreeTime](https://github.com/neherlab/treetime).
There are _many_ options that can be specified here in `refine` to help you get a good tree.
You can read about them in [the bioinformatics commands reference](/docs/bioinformatics/commands#refine).

`refine` will produce as output:
* another tree (newick format)
* a JSON format file with further information about each node

```
augur refine \
    --tree results/tree_raw.nwk \
    --alignment results/masked.vcf.gz \
    --vcf-reference data/ref.fasta \
    --metadata data/meta.tsv \
    --timetree \
    --root residual \
    --coalescent opt \
    --output-tree results/tree.nwk \
    --output-node-data results/branch_lengths.json
```

In addition to assigning times to internal nodes, the `refine` command filters tips that are likely outliers.
Branch lengths in the resulting Newick tree measure adjusted nucleotide divergence.
All other data inferred by TreeTime is stored by strain or internal node name in the JSON file.

## Annotate the Phylogeny

Now that we have an accurate tree and some information about the ancestral sequences, we can annotate some interesting data onto our phylogeny.
TreeTime can infer ancestral sequences and ancestral traits from an existing phylogenetic tree and metadata to annotate each tip of the tree.

### Infer Ancestral Sequences

We can reconstruct the ancestral sequences for the internal nodes on our phylogeny and identify any nucleotide mutations on the branches leading to any node in the tree.
You can read about all the options for `ancestral` in [the bioinformatics commands reference](/docs/bioinformatics/commands#ancestral).

For VCF runs, `ancestral` will produce another VCF that contains entries for the reconstructed sequence of all the internal nodes, as well as a JSON-format file that contains nucleotide mutation information for each node.

```
augur ancestral \
    --tree results/tree.nwk \
    --alignment results/masked.vcf.gz \
    --vcf-reference data/ref.fasta \
    --inference joint \
    --output results/nt_muts.json \
    --output-vcf results/nt_muts.vcf
```

### Identify Amino-Acid Mutations

With `translate` we can identify amino acid mutations from the nucleotide mutations and a GFF file with gene coordinate annotations.
The resulting JSON file contains amino acid mutations indexed by strain or internal node name and by gene name.
`translate` will also produce a VCF-style file with the amino acid changes for each gene and each sequence, and FASTA file with the translated 'reference' genes which the VCF-style file 'maps' to.

Because of the number of genes in TB, we will only translate some genes to save time. We can pass in a list of genes to translate (genes associated with drug resistance) using `--genes`.
Note that the `--reference-sequence` option is how you pass in the GFF file with the gene coordinates.

```
augur translate \
    --tree results/tree.nwk \
    --ancestral-sequences results/nt_muts.vcf \
    --vcf-reference data/ref.fasta \
    --genes config/genes.txt \
    --reference-sequence config/Mtb_H37Rv_NCBI_Annot.gff \
    --output results/aa_muts.json \
    --alignment-output results/translations.vcf \
    --vcf-reference-output results/translations_reference.fasta
```

### Reconstruct Ancestral States

`traits` can reconstruct the probable ancestral state of traits like location and host (or others).
This is done by specifying a column or columns in the meta-data file.

`--confidence` will give confidence estimates for the reconstructed states.
The output will be a JSON file with the state (and confidence, if specified) information for each node.

```
augur traits \
    --tree results/tree.nwk \
    --metadata data/meta.tsv \
    --columns location \
    --confidence \
    --output results/traits.json
```

### Identify Specified Clades

In the [original paper](http://www.pnas.org/content/112/44/13609), the authors identified 'sublineages' within the dataset.
We can add these to our dataset as 'clades' by defining the sublineages with amino-acid or nucleotide mutations specific to that sublineage, given here in the file `clades.tsv`.

As clades, these sublineages will be labelled and we'll be able to color the tree by them.

```
augur clades \
    --tree results/tree.nwk \
    --mutations results/nt_muts.json results/aa_muts.json \
    --clades config/clades.tsv \
    --output results/clades.json
```

### Identify Drug Resistance Mutations

`sequence-traits` can identify any trait associated with particular nucleotide or amino-acid mutations, not just drug resistance mutations.

This dataset doesn't actually contain any drug resistance mutations, but identifying such mutations is often of interest to those working on tuberculosis.
Here, we'll run this step as an example, even though it won't add anything to the tree for this dataset.

```
augur sequence-traits \
    --ancestral-sequences results/nt_muts.vcf \
    --vcf-reference data/ref.fasta \
    --translations results/translations.vcf \
    --vcf-translate-reference results/translations_reference.fasta \
    --features config/DRMs-AAnuc.tsv \
    --count traits \
    --label Drug_Resistance \
    --output results/drms.json
```

## Export the Results

Finally, collect all node annotations and metadata and export it all in auspice’s JSON format.
The resulting tree and metadata JSON files are the inputs to the auspice visualization tool.

```
augur export \
    --tree results/tree.nwk \
    --metadata data/meta.tsv \
    --node-data results/branch_lengths.json \
                results/traits.json \
                results/aa_muts.json \
                results/nt_muts.json \
                results/clades.json \
                results/drms.json \
    --auspice-config config/config.json \
    --colors config/color.tsv \
    --lat-longs config/lat_longs.tsv \
    --output-tree auspice/tb_tree.json \
    --output-meta auspice/tb_meta.json
```

## Visualize the Results

If you entered the Nextstrain build environment using `nextstrain shell` at the beginning of this tutorial, leave it now using the `exit` command and then use `nextstrain view` to visualize the TB build output in `auspice/*.json`.

```
# Leave the shell you entered earlier.
exit

# View results in your auspice/ directory.
nextstrain view auspice/
```

If you're not using the Nextstrain CLI shell, then copy the `auspice/*.json` files into the `data` directory of your local auspice installation and start auspice from there.
You can use the commands below (adjusted if necessary), or copy them using a graphical file explorer.

```
# Copy files into auspice data directory.  Adjust
# paths if auspice isn't installed in ~/src/auspice/.
mkdir ~/src/auspice/data/
cp auspice/*.json ~/src/auspice/data/
```

Then enter your `auspice` directry and start auspice.

```
# Navigate into auspice.
cd ~/src/auspice/data/

# Start auspice.
npm run dev
```

When auspice is running, navigate to <http://localhost:4000/local/tb> in your browser to view the results.

To stop auspice and return to the command line when you are done viewing your data, press CTRL+C.

## Automate the Build with Snakemake

While it is instructive to run all of the above commands manually, it is more practical to automate their execution with a single script.
Nextstrain implements these automated pathogen builds with [Snakemake](https://snakemake.readthedocs.io) by defining a `Snakefile` like [the one in the TB repository you downloaded](https://github.com/nextstrain/tb/blob/master/Snakefile).

First delete the output from the manual steps above.
(Be sure to navigate into the `tb-tutorial/` directory first.)

```
rm -rf results/ auspice/
```

Then, if you're using the Nextstrain CLI tool, run:

```
nextstrain build .
```

to run the automated pathogen build.

If you're not using the Nextstrain CLI tool, run:

```
snakemake
```

The automated build runs all of the manual steps above up through the auspice export.
View the results the same way you did before to confirm it produced the same TB build you made manually.

Note that automated builds will only re-run steps when the data changes.
This means builds will pick up where they left off if they are restarted after being interrupted.
If you want to force a re-run of the whole build, first remove any previous output with `nextstrain build . clean` or `snakemake clean`.

## Next steps

* Learn more about [augur commands](../bioinformatics).

* Learn more about [auspice visualizations](../visualisation).

* Learn more about [creating and modifying snakemake files](../pathogen-builds/snakemake).

* Fork the [TB pathogen repository on GitHub](https://github.com/nextstrain/tb), modify the Snakefile to make your own pathogen build, and learn [how to publish your build on nextstrain.org](../visualisation/introduction#viewing-your-data-through-nextstrainorg).
