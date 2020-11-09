---
title: "Augur: a bioinformatics toolkit for phylogenetic analyses of human pathogens"
papersize: a4
tags:
  - Python
  - bioinformatics
  - molecular biology
  - human pathogens
  - evolutionary biology
date: 5 November 2020
bibliography: paper.bib
---

# Summary and statement of need

The analysis of human pathogens requires a diverse collection of bioinformatics tools.
These tools include standard genomic and phylogenetic software and custom software developed to handle the relatively numerous and short genomes of viruses and bacteria.
Researchers increasingly depend on the outputs of these tools to infer transmission dynamics of human diseases and make actionable recommendations to public health officials [@Black2020; @Gardy2015].
Under these circumstances, bioinformatics tools must scale rapidly with the number of disease samples to enable real-time analyses of pathogen evolution.
To meet these needs, we developed Augur, a bioinformatics toolkit designed for phylogenetic analyses of human pathogens.

Augur originally existed as an internal component of the nextflu [@Neher:2015jr] and Nextstrain [@Hadfield2018] applications.
The original nextflu scripts only supported seasonal influenza viruses.
When nextflu was replaced with Nextstrain and expanded to support multiple viral and bacterial pathogens, each pathogen received its own copy of the original scripts.
The resulting redundancy of these large scripts complicated efforts to debug analyses, add new features for all pathogens, and add support for new pathogens.
In its original form, Augur consisted of two monolithic Python scripts, "prepare" and "process", that performed most operations in memory.
These scripts prepared a subset of pathogen sequences and metadata and then processed those data to produce an annotated phylogeny that could be viewed at [Nextstrain](https://nextstrain.org).
Critically, this software architecture led to long-lived, divergent branches of untested code in version control that Nextstrain team members could not confidently merge without potentially breaking existing analyses.

# Implementation

To address these issues, we refactored the original Augur scripts into a toolkit of individual subcommands wrapped by a single command line executable, `augur`.
With this approach, we followed the pattern established by samtools [@Li2009] and bcftools [@Li2011] where subcommands perform single, tightly-scoped tasks (e.g., "view", "sort", "merge", etc.) that can be chained together in bioinformatics pipelines.
We migrated or rewrote the existing functionality of the original Augur scripts into appropriate corresponding Augur subcommands.
To enable interoperability with existing bioinformatics tools, we designed subcommands to accept inputs and produce outputs in standard bioinformatics file formats wherever possible.
For example, we represented all raw sequence data in FASTA format, alignments in either FASTA or VCF format, and phylogenies in Newick format.
To handle the common case where a standard file format could not represent some or all of the outputs produced by an Augur command, we implemented a lightweight JSON schema to store the remaining data.
The "node data" JSON format represents one such Augur-specific file format that supports arbitrary annotations of phylogenies indexed by the name assigned to internal nodes or tips.
To provide a standard interface for our own analyses, we also designed several Augur subcommands to wrap existing bioinformatics tools including `augur align` (mafft [@Katoh2002]) and `augur tree` (FastTree [@Price2010], RAxML [@Stamatakis2014], and IQ-TREE [@Nguyen2014]).
Many commands including `augur refine`, `traits` and `ancestral` make extensive use of TreeTime [@Sagulenko2018] to provide time-scaled phylogenetic trees or further annotate the phylogeny.

By implementing the core components of Augur as a command line tool, we were able to rewrite our existing pathogen analyses as straightforward bioinformatics workflows using existing workflow management software like Snakemake [@Snakemake].
Most pathogen workflows begin with user-curated sequences in a FASTA file (e.g., `sequences.fasta`) and metadata describing each sequence in a tab-delimited text file (e.g., `metadata.tsv`).
Users can apply a series of Augur commands and other standard bioinformatics tools to these files to create annotated phylogenies that can be viewed in Auspice, the web application that serves [Nextstrain](https://nextstrain.org) (\autoref{fig:example-workflows}).
This approach allows users to leverage the distributed computing abilities of workflow managers to run multiple steps of the workflow in parallel and also run individual commands that support multiprocessing in parallel.

The modular Augur interface has enabled phylogenetic and genomic epidemiological analyses by academic researchers, public health laboratories, and private companies.
Most recently, these tools have supported the real-time tracking of SARS-CoV-2 evolution at global and local scales [@NextstrainNcov2020, @Bedford2020, @Alm2020].
This success has attracted contributions from the open source community that have allowed us to improve Augur's functionality, documentation, and test coverage.
Augur can be installed from PyPI ([nextstrain-augur](https://pypi.org/project/nextstrain-augur/)) and Bioconda ([augur](https://bioconda.github.io/recipes/augur/README.html)).
[See the full documentation](http://docs.nextstrain.org/) for more details about how to use or contribute to development of Augur.

# Figures

![Example workflows composed with Snakemake from Augur commands for A) Zika virus and B) tuberculosis.
  Each node in the workflow graph represents an Augur command than performs a specific part of the analysis (e.g., aligning sequences, building a tree, etc.).
  A typical workflow starts by filtering sequences and metadata to a desired subset for analysis followed by inference of a phylogeny, annotation of that phylogeny, and export of the annotated phylogeny to a JSON that can be viewed on Nextstrain.
  Workflows for viral (A) and bacterial (B) pathogens follow a similar structure but also support custom pathogen-specific steps.
  Multiple outgoing edges from a single node represent opportunities to run the workflow in parallel.
  See the full workflows at [https://github.com/nextstrain/zika-tutorial](https://github.com/nextstrain/zika-tutorial) and [https://github.com/nextstrain/tb](https://github.com/nextstrain/tb).\label{fig:example-workflows}](example-modular-augur-workflows.pdf)

# Acknowledgments

Thank you to all of [the open source community members who have contributed to Augur](https://github.com/nextstrain/augur/graphs/contributors).
Thank you to Dan Fornika from BCCDC Public Health Laboratory for creating the first conda recipe for Augur in Bioconda.

# References