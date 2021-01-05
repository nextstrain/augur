# This file is essentially a copy of ../zika/Snakefile
# but is here as it doesn't include the all / clean rule
# which, once imported, couldn't (?!?) be overridden.

rule base_files:
    params:
        input_fasta = "../zika/data/zika.fasta",
        dropped_strains = "../zika/config/dropped_strains.txt",
        reference = "../zika/config/zika_outgroup.gb",
        colors = "../zika/config/colors.tsv",
        auspice_config_v2 = "../zika/config/auspice_config_v2.json",

base_files = rules.base_files.params

rule parse:
    input:
        sequences = base_files.input_fasta
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/metadata.tsv"
    params:
        fasta_fields = "strain virus accession date region country division city db segment authors url title journal paper_url"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message: "Filtering sequences"
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        exclude = base_files.dropped_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country year month",
        sequences_per_group = 1,
        min_date = 2012
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --no-probabilistic-sampling \
            --min-date {params.min_date}
        """

rule align:
    message: "Aligning sequences"
    input:
        sequences = rules.filter.output.sequences,
        reference = base_files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
        """

rule refine:
    message: "Refining tree"
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """
