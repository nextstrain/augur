# What is a "build"?

Nextstrain's focus on providing a _real-time_ snapshot of evolving pathogen populations necessitates a reproducible analysis that can be rerun when new sequences are available.
The individual steps necessary to repeat analysis together comprise a "build".


Because no two datasets or pathogens are the same, we build Augur to be flexible and suitable for different analyses.
The individual Augur commands are composable, and can be mixed and matched with other scripts as needed.
These steps, taken together, are what we refer to as a "build".


### Example build

The [Zika virus tutorial](https://docs.nextstrain.org/en/latest/tutorials/creating-a-phylogenetic-workflow.html#run-a-nextstrain-build) describes a build which contains the following steps:

1. Prepare pathogen sequences and metadata
2. Align sequences
3. Construct a phylogeny from aligned sequences
4. Annotate the phylogeny with inferred ancestral pathogen dates, sequences, and traits
5. Export the annotated phylogeny and corresponding metadata into auspice-readable format

and each of these can be run via a separate `augur` command.

If you look at the [other tutorials](https://docs.nextstrain.org/en/latest/tutorials/index.html), each one uses a slightly different combination of `augur` commands depending on the pathogen.

### Snakemake

While it is possible to run a build by running each of the individual steps, we typically group these together into a make-type file.
[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) is "a tool to create reproducible and scalable data analyses... via a human-readable, Python-based language."

> Snakemake is installed as part of the [conda environment](https://docs.nextstrain.org/en/latest/guides/install/local-installation.html) or the [docker container](https://docs.nextstrain.org/en/latest/guides/install/cli-install.html).
If you ever see a build which has a "Snakefile" then you can run this by typing `snakemake --cores 1` or `nextstrain build --cpus 1 .`, respectively.
