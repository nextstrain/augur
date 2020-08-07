## Augur and snakemake

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
