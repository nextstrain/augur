import argparse

from .__version__ import __version__

command_strings = {
    "parse": "Parse delimited fields from FASTA sequence names into a TSV and FASTA file.",
    "filter": "Filter and subsample a sequence set.",
    "mask": "Mask specified sites from a VCF file.",
    "align": "Align multiple sequences from FASTA.",
    "tree": "Build a tree using a variety of methods.",
    "refine": "Refine an initial tree using sequence metadata.",
    "ancestral": "Infer ancestral sequences based on a tree.",
    "translate": "Translate gene regions from nucleotides to amino acids.",
    "reconstruct_sequences": "Reconstruct alignments from mutations inferred on the tree",
    "clades": "Assign clades to nodes in a tree based on amino-acid or nucleotide signatures.",
    "traits": "Infer ancestral traits based on a tree.",
    "sequence_traits": "Annotate sequences based on amino-acid or nucleotide signatures.",
    "lbi": "Calculate LBI for a given tree and one or more sets of parameters.",
    "distance": 'Calculate the distance between sequences across entire genes or at a predefined subset of sites.',
    "titers": "Annotate a tree with actual and inferred titer measurements.",
    "frequencies": "infer frequencies of mutations or clades",
    "export": "Export JSON files suitable for visualization with auspice.",
    "validate": "Validate files related to augur consumption or export.",
    "version": "Print the version of augur.",
    "import": "Import analyses into augur pipeline from other systems"
}

def augur_cli():
    parser = argparse.ArgumentParser(prog='augur',
                                     description="Augur: A bioinformatics toolkit for phylogenetic analysis.")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    _globals = globals()

    subparsers = parser.add_subparsers(dest='command')
    for command, help_str in command_strings.items():
        p = subparsers.add_parser(command, help=help_str, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        _globals[f'{command.replace("-", "_")}_options'](p)

    args = parser.parse_args()
    execute_run(args)

def execute_run(args):
    from importlib import import_module
    module = import_module(f'.{args.command}', package='augur')
    module.run(args)

def lbi_options(parser):
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--branch-lengths", help="JSON with branch lengths and internal node dates estimated by TreeTime", required=True)
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)
    parser.add_argument("--attribute-names", nargs="+", help="names to store distances associated with the corresponding masks", required=True)
    parser.add_argument("--tau", nargs="+", type=float, help="tau value(s) defining the neighborhood of each clade", required=True)
    parser.add_argument("--window", nargs="+", type=float, help="time window(s) to calculate LBI across", required=True)
    parser.add_argument("--no-normalization", action="store_true", help="disable normalization of LBI by the maximum value")

def frequencies_options(parser):
    parser.add_argument('--method', choices=["diffusion", "kde"], required=True,
                        help="method by which frequencies should be estimated")
    parser.add_argument('--metadata', type=str, required=True,
                        help="tab-delimited metadata including dates for given samples")
    parser.add_argument('--regions', type=str, nargs='+', default=['global'],
                        help="region to subsample to")
    parser.add_argument("--pivot-interval", type=int, default=3,
                        help="number of months between pivots")
    parser.add_argument('--min-date', type=float,
                        help="minimal pivot value")
    parser.add_argument('--max-date', type=float,
                        help="maximal pivot value")

    # Tree-specific arguments
    parser.add_argument('--tree', '-t', type=str,
                        help="tree to estimate clade frequencies for")
    parser.add_argument("--include-internal-nodes", action="store_true",
                        help="calculate frequencies for internal nodes as well as tips")

    # Alignment-specific arguments
    parser.add_argument('--alignments', type=str, nargs='+',
                        help="alignments to estimate mutations frequencies for")
    parser.add_argument('--gene-names', nargs='+', type=str,
                        help="names of the sequences in the alignment, same order assumed")
    parser.add_argument('--ignore-char', type=str, default='',
                        help="character to be ignored in frequency calculations")
    parser.add_argument('--minimal-frequency', type=float, default=0.05,
                        help="minimal all-time frequencies for a trajectory to be estimates")

    # KDE-specific arguments
    parser.add_argument("--narrow-bandwidth", type=float, default=1 / 12.0, help="the bandwidth for the narrow KDE")
    parser.add_argument("--wide-bandwidth", type=float, default=3 / 12.0, help="the bandwidth for the wide KDE")
    parser.add_argument("--proportion-wide", type=float, default=0.2,
                        help="the proportion of the wide bandwidth to use in the KDE mixture model")
    parser.add_argument("--weights",
                        help="a dictionary of key/value mappings in JSON format used to weight KDE tip frequencies")
    parser.add_argument("--weights-attribute",
                        help="name of the attribute on each tip whose values map to the given weights dictionary")
    parser.add_argument("--censored", action="store_true", help="calculate censored frequencies at each pivot")

    # Diffusion frequency specific arguments
    parser.add_argument('--minimal-clade-size', type=int, default=0,
                        help="minimal number of tips a clade must have for its diffusion frequencies to be reported")
    parser.add_argument('--minimal-clade-size-to-estimate', type=int, default=10,
                        help="""minimal number of tips a clade must have for its diffusion frequencies to be estimated
                                    by the diffusion likelihood; all smaller clades will inherit frequencies from their
                                    parents""")
    parser.add_argument("--stiffness", type=float, default=10.0,
                        help="parameter penalizing curvature of the frequency trajectory")
    parser.add_argument("--inertia", type=float, default=0.0, help="determines how frequencies continue "
                                                                   "in absense of data (inertia=0 -> go flat, inertia=1.0 -> continue current trend)")

    # Output arguments
    parser.add_argument('--output-format', default='auspice', choices=['auspice', 'nextflu'],
                        help="format to export frequencies JSON depending on the viewing interface")
    parser.add_argument('--output', '-o', type=str,
                        help='JSON file to save estimated frequencies to')


def export_options(parser):
    metavar_msg ="Augur export now needs you to define the JSON version " + \
                 "you want, e.g. `augur export v2`."
    subparsers = parser.add_subparsers(title="JSON SCHEMA",
                                       metavar=metavar_msg)
    subparsers.required = True

    v1 = subparsers.add_parser('v1', help="Export version 1 JSON schema (separate meta and tree JSONs)")
    core = v1.add_argument_group("REQUIRED")
    core.add_argument('--tree','-t', required=True, help="tree to perform trait reconstruction on")
    core.add_argument('--metadata', required=True, help="tsv file with sequence meta data")
    core.add_argument('--node-data', required=True, nargs='+', help="JSON files with meta data for each node")
    core.add_argument('--output-tree', help="JSON file name that is passed on to auspice (e.g., zika_tree.json).")
    core.add_argument('--output-meta', help="JSON file name that is passed on to auspice (e.g., zika_meta.json).")
    core.add_argument('--auspice-config', help="file with auspice configuration")

    options = v1.add_argument_group("OPTIONS")
    options.add_argument('--colors', help="file with color definitions")
    options.add_argument('--lat-longs', help="file latitudes and longitudes, overrides built in mappings")
    options.add_argument('--tree-name', default=False, help="Tree name (needed for tangle tree functionality)")
    options.add_argument('--minify-json', action="store_true", help="export JSONs without indentation or line returns")
    options.add_argument('--output-sequence', help="JSON file name that is passed on to auspice (e.g., zika_seq.json).")
    options.add_argument('--reference', required=False, help="reference sequence for export to browser, only vcf")
    options.add_argument('--reference-translations', required=False, help="reference translations for export to browser, only vcf")

    v1.add_argument("--v1", help=argparse.SUPPRESS, default=True)

    v2 = subparsers.add_parser("v2", help="Export version 2 JSON schema")

    required = v2.add_argument_group(
        title="REQUIRED"
    )
    required.add_argument('--tree', '-t', metavar="newick", required=True,
                          help="Phylogenetic tree, usually output from `augur refine`")
    required.add_argument('--node-data', metavar="JSON", required=True, nargs='+',
                          help="JSON files containing metadata for nodes in the tree")
    required.add_argument('--output', metavar="JSON", required=True,
                          help="Ouput file (typically for visualisation in auspice)")

    config = v2.add_argument_group(
        title="DISPLAY CONFIGURATION",
        description="These control the display settings for auspice. \
                You can supply a config JSON (which has all available options) or command line arguments (which are more limited but great to get started). \
                Supplying both is fine too, command line args will overrule what is set in the config file!"
    )
    config.add_argument('--auspice-config', metavar="JSON", help="Auspice configuration file")
    config.add_argument('--title', type=str, metavar="title", help="Title to be displayed by auspice")
    config.add_argument('--maintainers', metavar="name", action="append", nargs='+',
                        help="Analysis maintained by, in format 'Name <URL>' 'Name2 <URL>', ...")
    config.add_argument('--build-url', type=str, metavar="url", help="Build URL/repository to be displayed by Auspice")
    config.add_argument('--description', metavar="description.md",
                        help="Markdown file with description of build and/or acknowledgements to be displayed by Auspice")
    config.add_argument('--geo-resolutions', metavar="trait", nargs='+',
                        help="Geographic traits to be displayed on map")
    config.add_argument('--color-by-metadata', metavar="trait", nargs='+',
                        help="Metadata columns to include as coloring options")
    config.add_argument('--panels', metavar="panels", nargs='+', choices=['tree', 'map', 'entropy', 'frequencies'],
                        help="Restrict panel display in auspice. Options are %(choices)s. Ignore this option to display all available panels.")

    optional_inputs = v2.add_argument_group(
        title="OPTIONAL INPUT FILES"
    )
    optional_inputs.add_argument('--metadata', metavar="TSV", help="Additional metadata for strains in the tree")
    optional_inputs.add_argument('--colors', metavar="TSV", help="Custom color definitions")
    optional_inputs.add_argument('--lat-longs', metavar="TSV",
                                 help="Latitudes and longitudes for geography traits (overrides built in mappings)")

    optional_settings = v2.add_argument_group(
        title="OPTIONAL SETTINGS"
    )
    optional_settings.add_argument('--minify-json', action="store_true",
                                   help="export JSONs without indentation or line returns")
    optional_settings.add_argument('--include-root-sequence', action="store_true",
                                   help="Export an additional JSON containing the root sequence (reference sequence for vcf) used to identify mutations. The filename will follow the pattern of <OUTPUT>_root-sequence.json for a main auspice JSON of <OUTPUT>.json")


def validate_options(parser):
    subparsers = parser.add_subparsers(dest="subcommand", help="Which file(s) do you want to validate?")

    subparsers.add_parser("export-v2", help="validate JSON intended for auspice v2") \
        .add_argument('main_json', metavar='JSON', help="exported (main) v2 auspice JSON")

    export_v1 = subparsers.add_parser("export-v1", help="validate tree+meta JSONs intended for auspice v1")
    export_v1.add_argument('meta_json', metavar='META-JSON', help="exported (v1) meta JSON")
    export_v1.add_argument('tree_json', metavar='TREE-JSON', help="exported (v1) tree JSON")

    subparsers.add_parser("auspice-config-v2", help="validate auspice config intended for `augur export v2`") \
        .add_argument('config_json', metavar='JSON', help="auspice config JSON")

def version_options(parser):
    pass

def import_options(parser):
    metavar_msg = "Import analyses into augur pipeline from other systems"
    subparsers = parser.add_subparsers(title="TYPE",
                                       metavar=metavar_msg)
    subparsers.required = True

    beast_parser = subparsers.add_parser('beast', help="Import beast analysis")
    beast_parser.add_argument("--beast", help=argparse.SUPPRESS, default=True) # used to disambiguate subcommands
    beast_parser.add_argument('--mcc', required=True, help="BEAST MCC tree")
    beast_parser.add_argument('--most-recent-tip-date', default=0, type=float, help='Numeric date of most recent tip in tree (--tip-date-regex, --tip-date-format and --tip-date-delimeter are ignored if this is set)')
    beast_parser.add_argument('--tip-date-regex', default=r'[0-9]{4}(\-[0-9]{2})*(\-[0-9]{2})*$', type=str, help='regex to extract dates from tip names')
    beast_parser.add_argument('--tip-date-format', default="%Y-%m-%d", type=str, help='Format of date (if extracted by regex)')
    beast_parser.add_argument('--tip-date-delimeter', default="-", type=str, help='delimeter used in tip-date-format. Used to match partial dates.')
    beast_parser.add_argument('--verbose', action="store_true", help="Display verbose output. Only useful for debugging.")
    beast_parser.add_argument('--recursion-limit', default=False, type=int, help="Set a custom recursion limit (dangerous!)")
    beast_parser.add_argument('--output-tree', required=True, type=str, help='file name to write tree to')
    beast_parser.add_argument('--output-node-data', required=True, type=str, help='file name to write (temporal) branch lengths & BEAST traits as node data')

def distance_options(parser):
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str,
                        help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--attribute-name", nargs="+",
                        help="name to store distances associated with the given distance map; multiple attribute names are linked to corresponding positional comparison method and distance map arguments",
                        required=True)
    parser.add_argument("--compare-to", nargs="+", choices=["root", "ancestor", "pairwise"],
                        help="type of comparison between samples in the given tree including comparison of all nodes to the root (root), all tips to their last ancestor from a previous season (ancestor), or all tips from the current season to all tips in previous seasons (pairwise)",
                        required=True)
    parser.add_argument("--map", nargs="+",
                        help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates",
                        required=True)
    parser.add_argument("--date-annotations",
                        help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date")
    parser.add_argument("--earliest-date",
                        help="earliest date at which samples are considered to be from previous seasons (e.g., 2019-01-01). This date is only used in pairwise comparisons. If omitted, all samples prior to the latest date will be considered.")
    parser.add_argument("--latest-date",
                        help="latest date at which samples are considered to be from previous seasons (e.g., 2019-01-01); samples from any date after this are considered part of the current season")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name",
                        required=True)

def titers_options(parser):
    subparsers = parser.add_subparsers()

    tree_model = subparsers.add_parser('tree', help='tree model')
    tree_model.add_argument('--titers', nargs='+', type=str, required=True, help="file with titer measurements")
    tree_model.add_argument('--tree', '-t', type=str, required=True, help="tree to perform fit titer model to")
    tree_model.add_argument('--allow-empty-model', action="store_true", help="allow model to be empty")
    tree_model.add_argument('--output', '-o', type=str, required=True, help='JSON file to save titer model')

    sub_model = subparsers.add_parser('sub', help='substitution model')
    sub_model.add_argument('--titers', nargs='+', type=str, required=True, help="file with titer measurements")
    sub_model.add_argument('--alignment', nargs='+', type=str, required=True,
                           help="sequence to be used in the substitution model, supplied as fasta files")
    sub_model.add_argument('--gene-names', nargs='+', type=str, required=True,
                           help="names of the sequences in the alignment, same order assumed")
    sub_model.add_argument('--tree', '-t', type=str, help="optional tree to annotate fit titer model to")
    sub_model.add_argument('--allow-empty-model', action="store_true", help="allow model to be empty")
    sub_model.add_argument('--output', '-o', type=str, required=True, help='JSON file to save titer model')

def ancestral_options(parser):
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--output-node-data', type=str,
                        help='name of JSON file to save mutations and ancestral sequences to')
    parser.add_argument('--output', '-o', type=str, help='DEPRECATED. Same as --output-node-data')
    parser.add_argument('--output-sequences', type=str,
                        help='name of FASTA file to save ancestral sequences to (FASTA alignments only)')
    parser.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                        help="calculate joint or marginal maximum likelihood ancestral sequence states")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')
    ambiguous = parser.add_mutually_exclusive_group()
    ambiguous.add_argument('--keep-ambiguous', action="store_false", dest='infer_ambiguous',
                           help='do not infer nucleotides at ambiguous (N) sites on tip sequences (leave as N).')
    ambiguous.add_argument('--infer-ambiguous', action="store_true",
                           help='infer nucleotides at ambiguous (N,W,R,..) sites on tip sequences and replace with most likely state.')
    parser.add_argument('--keep-overhangs', action="store_true", default=False,
                        help='do not infer nucleotides for gaps (-) on either side of the alignment')

def translate_options(parser):
    parser.add_argument('--tree', help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--ancestral-sequences', type=str,
                        help='JSON (fasta input) or VCF (VCF input) containing ancestral and tip sequences')
    parser.add_argument('--reference-sequence', required=True,
                        help='GenBank or GFF file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate (list or file containing list)")
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save aa-mutations to')
    parser.add_argument('--output', '-o', type=str, help='DEPRECATED. Same as --output-node-data')
    parser.add_argument('--alignment-output', type=str, help="write out translated gene alignments. "
                                                             "If a VCF-input, a .vcf or .vcf.gz will be output here (depending on file ending). If fasta-input, specify the file name "
                                                             "like so: 'my_alignment_%%GENE.fasta', where '%%GENE' will be replaced by the name of the gene")
    parser.add_argument('--vcf-reference-output', type=str,
                        help="fasta file where reference sequence translations for VCF input will be written")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')


def reconstruct_sequences_options(parser):
    parser.add_argument('--tree', required=True, help="tree as Newick file")
    parser.add_argument('--gene', type=str, help="gene to translate (list or file containing list)")
    parser.add_argument('--mutations', required=True, type=str, help="json file containing mutations "
                                                                     "mapped to each branch and the sequence of the root.")
    parser.add_argument('--vcf-aa-reference', type=str,
                        help='fasta file of the reference gene translations for VCF format')
    parser.add_argument('--internal-nodes', action='store_true', help="include sequences of internal nodes in output")
    parser.add_argument('--output', type=str)


def clades_options(parser):
    parser.add_argument('--tree', help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--mutations', nargs='+',
                        help='JSON(s) containing ancestral and tip nucleotide and/or amino-acid mutations ')
    parser.add_argument('--reference', nargs='+',
                        help='fasta files containing reference and tip nucleotide and/or amino-acid sequences ')
    parser.add_argument('--clades', type=str, help='TSV file containing clade definitions by amino-acid')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save clade assignments to')
    parser.add_argument('--output', '-o', type=str, help='DEPRECATED. Same as --output-node-data')

def traits_options(parser):
    parser.add_argument('--tree', '-t', required=True, help="tree to perform trait reconstruction on")
    parser.add_argument('--metadata', required=True, help="tsv/csv table with meta data")
    parser.add_argument('--weights', required=False,
                        help="tsv/csv table with equilibrium probabilities of discrete states")
    parser.add_argument('--columns', required=True, nargs='+',
                        help='metadata fields to perform discrete reconstruction on')
    parser.add_argument('--confidence', action="store_true",
                        help='record the distribution of subleading mugration states')
    parser.add_argument('--sampling-bias-correction', type=float,
                        help='a rough estimate of how many more events would have been observed'
                             ' if sequences represented an even sample. This should be'
                             ' roughly the (1-sum_i p_i^2)/(1-sum_i t_i^2), where p_i'
                             ' are the equilibrium frequencies and t_i are apparent ones.'
                             '(or rather the time spent in a particular state on the tree)')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save trait inferences to')
    parser.add_argument('--output', '-o', type=str, help='DEPRECATED. Same as --output-node-data')
    parser.epilog = "Note that missing data must be represented by a `?` character. Missing data will currently be inferred."

def sequence_traits_options(parser):
    parser.add_argument('--ancestral-sequences', type=str,
                        help="nucleotide alignment (VCF) to search for sequence traits in (can be generated from 'ancestral' using '--output-vcf')")
    parser.add_argument('--translations', type=str,
                        help="AA alignment to search for sequence traits in (can include ancestral sequences)")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the nucleotide VCF was mapped to')
    parser.add_argument('--vcf-translate-reference', type=str,
                        help='fasta file of the sequence the translated VCF was mapped to')
    parser.add_argument('--features', type=str,
                        help='file that specifies sites defining the features in a tab-delimited format: "GENE SITE ALT DISPLAY_NAME FEATURE". For nucleotide sites, GENE can be "nuc" (or column excluded entirely for all-nuc sites). "DISPLAY_NAME" can be blank or excluded entirely.')
    parser.add_argument('--count', type=str, choices=['traits', 'mutations'], default='traits',
                        help='Whether to count traits (ex: # drugs resistant to) or mutations')
    parser.add_argument('--label', type=str, default="# Traits", help='How to label the counts (ex: Drug_Resistance)')
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save sequence features to')
    parser.add_argument('--output', '-o', type=str, help='DEPRECATED. Same as --output-node-data')

def refine_options(parser):
    parser.add_argument('--alignment', '-a', help="alignment in fasta or VCF format")
    parser.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    parser.add_argument('--metadata', type=str, help="tsv/csv table with meta data for sequences")
    parser.add_argument('--output-tree', type=str, help='file name to write tree to')
    parser.add_argument('--output-node-data', type=str, help='file name to write branch lengths as node data')
    parser.add_argument('--timetree', action="store_true", help="produce timetree using treetime")
    parser.add_argument('--coalescent',
                        help="coalescent time scale in units of inverse clock rate (float), optimize as scalar ('opt'), or skyline ('skyline')")
    parser.add_argument('--gen-per-year', default=50, type=float,
                        help="number of generations per year, relevant for skyline output('skyline')")
    parser.add_argument('--clock-rate', type=float, help="fixed clock rate")
    parser.add_argument('--clock-std-dev', type=float, help="standard deviation of the fixed clock_rate estimate")
    parser.add_argument('--root', nargs="+", default='best',
                        help="rooting mechanism ('best', least-squares', 'min_dev', 'oldest') "
                             "OR node to root by OR two nodes indicating a monophyletic group to root by. "
                             "Run treetime -h for definitions of rooting methods.")
    parser.add_argument('--keep-root', action="store_true", help="do not reroot the tree; use it as-is. "
                                                                 "Overrides anything specified by --root.")
    parser.add_argument('--covariance', dest='covariance', action='store_true',
                        help="Account for covariation when estimating "
                             "rates and/or rerooting. "
                             "Use --no-covariance to turn off.")
    parser.add_argument('--no-covariance', dest='covariance',
                        action='store_false')  # If you set help here, it displays 'default: True' - which is confusing!
    parser.add_argument('--keep-polytomies', action='store_true', help='Do not attempt to resolve polytomies')
    parser.add_argument('--date-format', default="%Y-%m-%d", help="date format")
    parser.add_argument('--date-confidence', action="store_true", help="calculate confidence intervals for node dates")
    parser.add_argument('--date-inference', default='joint', choices=["joint", "marginal"],
                        help="assign internal nodes to their marginally most likely dates, not jointly most likely")
    parser.add_argument('--branch-length-inference', default='auto', choices=['auto', 'joint', 'marginal', 'input'],
                        help='branch length mode of treetime to use')
    parser.add_argument('--clock-filter-iqd', type=float, help='clock-filter: remove tips that deviate more than n_iqd '
                                                               'interquartile ranges from the root-to-tip vs time regression')
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--year-bounds', type=int, nargs='+',
                        help='specify min or max & min prediction bounds for samples with XX in year')
    parser.add_argument('--divergence-units', type=str, choices=['mutations', 'mutations-per-site'],
                        default='mutations-per-site', help='Units in which sequence divergences is exported.')
    parser.set_defaults(covariance=True)

def mask_options(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in VCF format")
    parser.add_argument('--mask', required=True, help="locations to be masked in BED file format")
    parser.add_argument('--output', '-o', help="output file")

def tree_options(parser):
    parser.add_argument('--alignment', '-a', required=True, help="alignment in fasta or VCF format")
    parser.add_argument('--method', default='iqtree', choices=["fasttree", "raxml", "iqtree"],
                        help="tree builder to use")
    parser.add_argument('--output', '-o', type=str, help='file name to write tree to')
    parser.add_argument('--substitution-model', default="GTR", choices=["HKY", "GTR", "HKY+G", "GTR+G", "GTR+R10"],
                        help='substitution model to use. Specify \'none\' to run ModelTest. Currently, only available for IQTREE.')
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                        help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    parser.add_argument('--exclude-sites', type=str,
                        help='file name of one-based sites to exclude for raw tree building (BED format in .bed files, DRM format in tab-delimited files, or one position per line)')
    parser.add_argument('--tree-builder-args', type=str, default='',
                        help='extra arguments to be passed directly to the executable of the requested tree method (e.g., --tree-builder-args="-czb")')


def filter_options(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--metadata', required=True, help="metadata associated with sequences")
    parser.add_argument('--min-date', type=float, help="minimal cutoff for numerical date")
    parser.add_argument('--max-date', type=float, help="maximal cutoff for numerical date")
    parser.add_argument('--min-length', type=int, help="minimal length of the sequences")
    parser.add_argument('--non-nucleotide', action='store_true',
                        help="exclude sequences that contain illegal characters")
    parser.add_argument('--exclude', type=str, help="file with list of strains that are to be excluded")
    parser.add_argument('--include', type=str,
                        help="file with list of strains that are to be included regardless of priorities or subsampling")
    parser.add_argument('--priority', type=str, help="file with list priority scores for sequences (strain\tpriority)")
    parser.add_argument('--sequences-per-group', type=int,
                        help="subsample to no more than this number of sequences per category")
    parser.add_argument('--group-by', nargs='+',
                        help="categories with respect to subsample; two virtual fields, \"month\" and \"year\", are supported if they don't already exist as real fields but a \"date\" field does exist")
    parser.add_argument('--subsample-seed',
                        help="random number generator seed to allow reproducible sub-sampling (with same input data). Can be number or string.")
    parser.add_argument('--exclude-where', nargs='+',
                        help="Exclude samples matching these conditions. Ex: \"host=rat\" or \"host!=rat\". Multiple values are processed as OR (matching any of those specified will be excluded), not AND")
    parser.add_argument('--include-where', nargs='+',
                        help="Include samples with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. This rule is applied last and ensures any sequences matching these rules will be included.")
    parser.add_argument('--output', '-o', help="output file", required=True)

def parse_options(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--output-sequences', help="output sequences file")
    parser.add_argument('--output-metadata', help="output metadata file")
    parser.add_argument('--fields', nargs='+', help="fields in fasta header")
    parser.add_argument('--prettify-fields', nargs='+',
                        help="apply string prettifying operations (underscores to spaces, capitalization, etc) to specified metadata fields")
    parser.add_argument('--separator', default='|', help="separator of fasta header")
    parser.add_argument('--fix-dates', choices=['dayfirst', 'monthfirst'],
                        help="attempt to parse non-standard dates and output them in standard YYYY-MM-DD format")

def align_options(parser):
    parser.add_argument('sequences', nargs="+", metavar="FASTA", help="sequences to align")
    parser.add_argument('--output', '-o', default="alignment.fasta", help="output file (default: %(default)s)")
    parser.add_argument('--nthreads', type=nthreads_value, default=1,
                        help="number of threads to use; specifying the value 'auto' will cause the number of available CPU cores on your system, if determinable, to be used")
    parser.add_argument('--method', default='mafft', choices=["mafft"], help="alignment program to use")
    parser.add_argument('--reference-name', metavar="NAME", type=str,
                        help="strip insertions relative to reference sequence; use if the reference is already in the input sequences")
    parser.add_argument('--reference-sequence', metavar="PATH", type=str,
                        help="Add this reference sequence to the dataset & strip insertions relative to this. Use if the reference is NOT already in the input sequences")
    parser.add_argument('--remove-reference', action="store_true", default=False,
                        help="remove reference sequence from the alignment")
    parser.add_argument('--fill-gaps', action="store_true", default=False,
                        help="If gaps represent missing data rather than true indels, replace by N after aligning.")
    parser.add_argument('--existing-alignment', metavar="FASTA", default=False,
                        help="An existing alignment to which the sequences will be added. The ouput alignment will be the same length as this existing alignment.")
    parser.add_argument('--debug', action="store_true", default=False,
                        help="Produce extra files (e.g. pre- and post-aligner files) which can help with debugging poor alignments.")

def available_cpu_cores(fallback: int=1) -> int:
    """
    Returns the number (an int) of CPU cores available to this **process**, if
    determinable, otherwise the number of CPU cores available to the
    **computer**, if determinable, otherwise the *fallback* number (which
    defaults to 1).
    """
    import os
    try:
        # Note that this is the correct function to use, not os.cpu_count(), as
        # described in the latter's documentation.
        #
        # The reason, which the documentation does not detail, is that
        # processes may be pinned or restricted to certain CPUs by setting
        # their "affinity".  This is not typical except in high-performance
        # computing environments, but if it is done, then a computer with say
        # 24 total cores may only allow our process to use 12.  If we tried to
        # naively use all 24, we'd end up with two threads across the 12 cores.
        # This would degrade performance rather than improve it!
        return len(os.sched_getaffinity(0))
    except:
        # cpu_count() returns None if the value is indeterminable.
        return os.cpu_count() or fallback


def nthreads_value(value):
    """
    Argument value validation and casting function for --nthreads.
    """

    if value.lower() == 'auto':
        return available_cpu_cores()

    try:
        return int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("'%s' is not an integer or the word 'auto'" % value) from None

if __name__ == "__main__":
    augur_cli()