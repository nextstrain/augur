"""
Infer ancestral sequences based on a tree.

The ancestral sequences are inferred using `TreeTime <https://academic.oup.com/ve/article/4/1/vex042/4794731>`_.
Each internal node gets assigned a nucleotide sequence that maximizes a
likelihood on the tree given its descendants and its parent node.
Each node then gets assigned a list of nucleotide mutations for any position
that has a mismatch between its own sequence and its parent's sequence.
The node sequences and mutations are output to a node-data JSON file.

If amino acid options are provided, the ancestral amino acid sequences for each
requested gene are inferred with the same method as the nucleotide sequences described above.
The inferred amino acid mutations will be included in the output node-data JSON
file, with the format equivalent to the output of `augur translate`.

The nucleotide and amino acid sequences are inferred separately in this command,
which can potentially result in mismatches between the nucleotide and amino
acid mutations. If you want amino acid mutations based on the inferred
nucleotide sequences, please use `augur translate`.

.. note::

    The mutation positions in the node-data JSON are one-based.
"""
from augur.argparse_ import ExtendOverwriteDefault
import argparse
from augur.errors import AugurError
import sys
import numpy as np
from Bio import SeqIO
from Bio.Phylo.BaseTree import Tree  # type: ignore[import-untyped]
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .translate import safe_translate
from .utils import parse_genes_argument, read_tree, InvalidTreeError, write_augur_json, get_json_name, \
    genome_features_to_auspice_annotation
from .io.file import open_file
from .io.sequences import read_single_sequence, is_vcf as is_filename_vcf
from treetime.vcf_utils import read_vcf, write_vcf
from collections import defaultdict
from .argparse_ import add_validation_arguments
from .util_support.node_data_file import NodeDataObject
from typing import cast, Any, Callable, TypedDict
from typing_extensions import NotRequired
from math import floor
from dataclasses import dataclass

VCF_Alignment = dict[str, dict[int, str]]

VCF_Metadata = dict[str, Any]

class Mutations(TypedDict):
    nodes: Any
    mask: Any # numpy.ndarray(bool)
    
class Ancestral_Reconstruction(TypedDict):
    tt: Any
    root_seq: str
    mutations: Mutations

class Nuc_Annotation(TypedDict):
    start: int
    end: int
    strand: str
    type: str

class Annotations_JSON(TypedDict):
    nuc: NotRequired[Nuc_Annotation]

class Ancestral_JSON(TypedDict):
    reference: dict[str, str]
    mask: NotRequired[str]
    annotations: Annotations_JSON
    nodes: Any

GENE_PATTERN = "%GENE"
"""String pattern used for gene replacement in filenames etc"""

def _make_seq_corrector(alphabet: str) -> Callable[[str], str]:
    """Build a function that replaces invalid or fully-ambiguous characters in a
    sequence with the standard ambiguous character for the given alphabet
    ('N' for nuc, 'X' for aa).

    A character is replaced if it is either:
    - not present in TreeTime's profile_map (i.e. completely unknown), or
    - fully ambiguous (its profile has equal weight across all states)

    The standard ambiguous character itself (N/X) is kept as-is since it is
    valid and expected by downstream code such as create_mask.
    
    The corrected sequence string will be all uppercase.
    """
    from treetime.seq_utils import alphabets, profile_maps
    profile_map = profile_maps[alphabet]
    n_states = len(alphabets[alphabet])

    if alphabet == 'nuc':
        ambiguous_char = 'N'
    elif alphabet == 'aa':
        ambiguous_char = 'X'
    else:
        raise ValueError(f"Unknown alphabet: {alphabet!r}")

    # Identify characters in the profile_map that are valid and not fully
    # ambiguous (partially ambiguous IUPAC chars like R, Y are kept)
    valid = set()
    for char, profile in profile_map.items():
        if char == ambiguous_char:
            valid.add(char)
        elif sum(profile) < n_states:
            valid.add(char)
        # else: fully ambiguous non-standard char → will be replaced

    def correct(seq: str) -> str:
        seq_upper = seq.upper()
        if all(c in valid for c in seq_upper):
            return seq_upper
        return ''.join(c if c in valid else ambiguous_char for c in seq_upper)

    return correct


def correct_alignment(aln_fname: str, correct_seq: Callable[[str], str]) -> MultipleSeqAlignment:
    """Read an alignment from a FASTA file and correct sequences using the
    provided correction function (from _make_seq_corrector).

    Returns a MultipleSeqAlignment suitable for passing directly to TreeAnc.
    """
    from Bio import AlignIO
    alignment = AlignIO.read(aln_fname, 'fasta')
    corrected_records = [
        SeqRecord(Seq(correct_seq(str(record.seq))), id=record.id, description=record.description)
        for record in alignment
    ]
    return MultipleSeqAlignment(corrected_records)


def create_mask(is_vcf, tt, reference_sequence, aln):
    """
    Identify sites for which every terminal sequence is ambiguous. These sites
    will be masked to prevent rounding errors in the maximum likelihood
    inference from assigning an arbitrary nucleotide to sites at internal nodes.

    Parameters
    ----------
    is_vcf : bool
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction. Unused if is_vcf.
    reference_sequence : str
        only used if is_vcf
    aln : dict
        describes variation (relative to reference) per sample. Only used if is_vcf.

    Returns
    -------
    numpy.ndarray(bool)
    """
    num_tips = len(tt.tree.get_terminals())
    ambiguous_count = np.zeros(tt.sequence_length, dtype=int)
    if is_vcf:
        variable_sites = set()
        # VCF ambiguous positions come in two forms:
        # Firstly, if VCF defines a "N" ALT and assigns it to every sample:
        for sample_data in aln.values():
            ambig_positions = []
            for pos, alt in sample_data.items():
                variable_sites.add(pos)
                if alt=='N':
                    ambig_positions.append(pos)
            # ambig_positions = [pos for pos, alt in sample_data.items() if alt=='N']
            np.add.at(ambiguous_count, ambig_positions, 1)
        # Secondly, if the VCF defines no mutations but the ref is "N":
        no_info_sites = np.array(list(reference_sequence)) == 'N'
        no_info_sites[list(variable_sites)] = False
        # and the mask is the union of these two forms
        mask = ambiguous_count==num_tips
        mask[no_info_sites] = True
    else:
        for n in tt.tree.get_terminals():
            ambiguous_count += np.array(tt.sequence(n,reconstructed=False, as_string=False)==tt.gtr.ambiguous, dtype=int)
        mask = ambiguous_count==num_tips
    return mask

def collect_mutations(tt, mask, reference_sequence=None, infer_ambiguous=False):
    """iterates of the tree and produces dictionaries with
    mutations and sequences for each node.

    If a reference sequence is provided then mutations can be collected for the
    root node. Masked positions at the root-node will be treated specially: if
    we infer ambiguity, then we report no mutations (i.e. we assume the
    reference base holds), otherwise we'll report a mutation from the <ref> to
    "N".

    Parameters
    ----------
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction
    mask : numpy.ndarray(bool)
    reference_sequence : str, optional

    Returns
    -------
    dict
        dict -> <node_name> -> [mut, mut, ...] where mut is a string in the form
        <from><1-based-pos><to>
    """

    data = {}
    inc = 1 # convert python numbering to start-at-1

    # Masked positions are fully ambiguous columns where the inferred base is
    # arbitrary; skip any mutations at those positions. This matters when
    # infer_ambiguous=False because tips retain their ambiguous character while
    # internal nodes carry an inferred base, producing spurious <inferred>→N
    # transitions in n.mutations.
    for n in tt.tree.find_clades():
        data[n.name] = [a+str(int(pos)+inc)+d
                        for a,pos,d in n.mutations if not mask[int(pos)]]

    if reference_sequence:
        data[tt.tree.root.name] = []
        for pos, (root_state, tree_state) in enumerate(zip(reference_sequence, tt.sequence(tt.tree.root, reconstructed=infer_ambiguous, as_string=True))):
            if mask[pos]:
                continue
            if root_state != tree_state:
                data[tt.tree.root.name].append(f"{root_state}{pos+1}{tree_state}")

    return data

def collect_sequences(tt, mask, reference_sequence=None, infer_ambiguous=False):
    """
    Create a full sequence for every node on the tree. Masked positions will
    have the reference base if we are inferring ambiguity, or the ambiguous
    character 'N'.

    Parameters
    ----------
    tt : treetime.TreeTime
        instance of treetime with valid ancestral reconstruction
    mask : numpy.ndarray(bool)
        Mask these positions by changing them to the ambiguous nucleotide
    reference_sequence : str or None
    infer_ambiguous : bool, optional
        if true, request the reconstructed sequences from treetime, otherwise retain input ambiguities

    Returns
    -------
    dict
        dict -> <node_name> -> sequence_string
    """
    sequences = {}

    ref_mask = None

    if reference_sequence and infer_ambiguous:
        ref_mask = np.array(list(reference_sequence))[mask]

    for n in tt.tree.find_clades():
        try:
            tmp = tt.sequence(n,reconstructed=infer_ambiguous, as_string=False)
            if ref_mask is not None:
                tmp[mask] = ref_mask
            else:
                tmp[mask] = tt.gtr.ambiguous
            sequences[n.name] = "".join(tmp)
        except:
            print("No sequence available for node ",n.name)
    return sequences


def run_ancestral(
        T:Tree,
        aln: str|VCF_Alignment|MultipleSeqAlignment,
        reference_sequence:str|None=None,
        is_vcf=False,
        full_sequences=False,
        fill_overhangs=False,
        infer_ambiguous=False,
        marginal=False,
        alphabet='nuc',
        rng_seed:int|None=None
        ) -> Ancestral_Reconstruction:
    """
    ancestral nucleotide reconstruction using TreeTime
    """
    from treetime import TreeAnc
    
    tt = TreeAnc(tree=T, aln=aln, ref=reference_sequence if is_vcf else None, gtr='JC69', alphabet=alphabet,
                 fill_overhangs=fill_overhangs, verbose=1, rng_seed=rng_seed)

    # only infer ancestral sequences, leave branch length untouched
    tt.infer_ancestral_sequences(infer_gtr=True, marginal=marginal,
                                 reconstruct_tip_states=infer_ambiguous, sample_from_profile='root')

    # add reference sequence to json structure. This is the sequence with
    # respect to which mutations on the tree are defined.
    if reference_sequence:
        root_seq = reference_sequence
    else:
        root_seq = str(tt.sequence(T.root, as_string=True))

    mask = create_mask(is_vcf, tt, reference_sequence, aln)
    mutations = collect_mutations(tt, mask, reference_sequence, infer_ambiguous)
    sequences = {}
    if full_sequences:
        sequences = collect_sequences(tt, mask, reference_sequence, infer_ambiguous)

    # Combine the mutations & sequences into a single dict which downstream code
    # expects
    nodes: defaultdict[str, dict[str, Any]] = defaultdict(dict)
    for n in tt.tree.find_clades():
        name = n.name
        if name in mutations:
            nodes[name]['muts'] = mutations[name]
        if name in sequences:
            nodes[name]['sequence'] = sequences[name]

    return {'tt': tt,
            'root_seq': root_seq,
            'mutations': {"nodes": nodes, "mask": mask}}


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("ancestral", help=__doc__)

    # ----------------------------- INPUTS ----------------------------- 
    input_group = parser.add_argument_group(
        "inputs",
        "Tree and sequences to use for ancestral reconstruction"
    )
    input_group.add_argument('--tree', '-t', required=True, help="prebuilt Newick")
    
    # ------------------- NUCLEOTIDE (TRANSLATION) OPTIONS ONLY ------------------------ 
    nucleotide_options_group = parser.add_argument_group(
        "nucleotide options",
        "Options to configure reconstruction of nucleotide sequences."
    )
    nucleotide_options_group.add_argument('--alignment', '-a', help="alignment in FASTA or VCF format")
    nucleotide_options_exclusive_group = nucleotide_options_group.add_mutually_exclusive_group()
    nucleotide_options_exclusive_group.add_argument('--vcf-reference', type=str, metavar='FASTA',
                                 help='[VCF alignment only] file of the sequence the VCF was mapped to.'
                                      ' Differences between this sequence and the inferred root will be reported as mutations on the root branch.')
    nucleotide_options_exclusive_group.add_argument('--root-sequence', type=str,metavar='FASTA/GenBank',
                                 help='[FASTA alignment only] file of the sequence that is used as root for mutation calling.'
                                      ' Differences between this sequence and the inferred root will be reported as mutations on the root branch.')

    # ----------------------------- GLOBAL OPTIONS -----------------------------
    global_options_group = parser.add_argument_group(
        "global options",
        "Options to configure reconstruction of both nucleotide and amino acid sequences"
    )
    global_options_group.add_argument('--inference', default='joint', choices=["joint", "marginal"],
                                    help="calculate joint or marginal maximum likelihood ancestral sequence states")
    global_options_group.add_argument('--seed', type=int, help="seed for random number generation")
    ambiguous = global_options_group.add_mutually_exclusive_group()
    ambiguous.add_argument('--keep-ambiguous', action="store_true",
                                help='do not infer ambiguous states on tip sequences')
    ambiguous.add_argument('--infer-ambiguous', action="store_true", default=True,
                                help='infer ambiguous states on tip sequences and replace with most likely state')
    global_options_group.add_argument('--keep-overhangs', action="store_true", default=False,
                                help='do not infer states for gaps (-) on either side of the alignment')

    # ------------------- AMINO-ACID (TRANSLATION) OPTIONS ONLY ------------------------ 
    amino_acid_options_group = parser.add_argument_group(
        "amino acid options",
        "Options to configure reconstruction of ancestral amino acid sequences."
    )
    amino_acid_options_group.add_argument('--genes', nargs='+', action=ExtendOverwriteDefault,
        help="gene(s) to translate (list or file containing list).")
    amino_acid_options_group.add_argument('--annotation',
                        help='GenBank or GFF file containing the annotation.')
    amino_acid_options_group.add_argument('--translations', type=str, help="Translated alignments for each CDS/Gene."
                           " If you are translating multiple genes you must specify the file name via a template"
                           f" like 'aa_sequences_%{GENE_PATTERN}.fasta' where %{GENE_PATTERN} will be replaced,"
                           " If you are translating a single gene using a pattern is optional."
                           " Currently only supported for FASTA-input.")
    amino_acid_options_group.add_argument('--report-inconsistent-translation', action="store_true",
                                help="Report where amino acid reconstruction differed from a translation of the reconstructed nuc sequence."
                                " Requires nucleotide reconstruction.")

    # ----------------------------- OUTPUTS ----------------------------- 
    output_group = parser.add_argument_group(
        "outputs",
        "Outputs supported for reconstructed ancestral sequences"
    )
    output_group.add_argument('--output-node-data', type=str, help='name of JSON file to save mutations and ancestral sequences to')
    output_group.add_argument('--output-sequences', type=str, help='name of FASTA file to save ancestral nucleotide sequences to (FASTA alignments only)')
    output_group.add_argument('--output-translations', type=str, help="name of the FASTA file(s) to save ancestral amino acid sequences to. "
                        f"Specify the file name via a template like 'ancestral_aa_sequences_%{GENE_PATTERN}.fasta' where %{GENE_PATTERN} will be replaced by"
                        "the gene name.")
    output_group.add_argument('--output-vcf', type=str, help='name of output VCF file which will include ancestral seqs')

    # ----------------------------- GENERAL / MISC ----------------------------- 
    general_group = parser.add_argument_group("general",)
    add_validation_arguments(general_group)

    return parser

@dataclass()
class Mode:
    is_vcf: bool
    nuc_reconstruction: bool
    aa_reconstruction: bool = False

def validate_arguments(args: argparse.Namespace) -> Mode:
    """
    Check that provided arguments are compatible.
    Where possible we use argparse built-ins, but they don't cover everything we want to check.
    This checking shouldn't be used by downstream code to assume arguments exist, however by checking for
    invalid combinations up-front we can exit quickly.
    """
    mode = Mode(
        nuc_reconstruction = bool(args.alignment),
        is_vcf = is_filename_vcf(args.alignment)
    )

    if not mode.nuc_reconstruction and args.root_sequence:
        raise AugurError("A root sequence is only used if you also reconstruct nuc sequences")

    # translations are run in one of two modes: multiple genes or a single gene
    required_aa_arguments = (args.annotation, args.genes, args.translations)
    if any(required_aa_arguments) and not all(required_aa_arguments):
        raise AugurError("For amino acid sequence reconstruction, you must provide an annotation file, a list of genes, and a path to amino acid sequences.")
    if all(required_aa_arguments):
        mode.aa_reconstruction = True

    if not mode.nuc_reconstruction and not mode.aa_reconstruction:
        raise AugurError("Neither nucleotide nor AA reconstruction requested")

    if args.report_inconsistent_translation and (not mode.aa_reconstruction or not mode.nuc_reconstruction):
        raise AugurError("Inconsistent translation reporting requires both nuc and AA reconstructions to be done")

    if not mode.nuc_reconstruction and args.output_sequences:
        raise AugurError("You cannot output (nuc) sequences if a (nuc) alignment is not provided")

    if args.output_sequences and args.output_vcf:
        raise AugurError("Both sequence (fasta) and VCF output have been requested, but these are incompatible.")

    if mode.is_vcf and args.output_sequences:
        raise AugurError("Sequence (fasta) output has been requested but the input alignment is VCF.")

    if not mode.is_vcf and args.output_vcf:
        raise AugurError("VCF output has been requested but the input alignment is not VCF.")

    return mode


def _read_tree(fname: str) -> Tree:
    try:
        T = read_tree(fname)
    except FileNotFoundError:
        raise AugurError(f"The provided tree file {fname!r} doesn't exist")
    except InvalidTreeError as error:
        raise AugurError(error)
    # Note that a number of other errors may be thrown by `read_tree` such as Bio.Phylo.NewickIO.NewickError

    missing_internal_node_names = [n.name is None for n in T.get_nonterminals()]
    if np.all(missing_internal_node_names):
        print("\n*** WARNING: Tree has no internal node names!", file=sys.stderr)
        print("*** Without internal node names, ancestral sequences can't be linked up to the correct node later.", file=sys.stderr)
        print("*** If you want to use 'augur export' or `augur translate` later, re-run this command with the output of 'augur refine'.", file=sys.stderr)
        print("*** If you haven't run 'augur refine', you can add node names to your tree by running:", file=sys.stderr)
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(fname) , file=sys.stderr)
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'", file=sys.stderr)
    
    return T

def _read_sequence_data(args: argparse.Namespace, is_vcf: bool) \
        -> tuple[str|None, VCF_Alignment|MultipleSeqAlignment, VCF_Metadata|None]:
    correct_nuc = _make_seq_corrector('nuc')
    if is_vcf:
        if not args.vcf_reference:
            raise AugurError("a reference Fasta is required with VCF-format alignments")
        compress_seq = read_vcf(args.alignment, args.vcf_reference)
        aln: VCF_Alignment | MultipleSeqAlignment = cast(VCF_Alignment, compress_seq['sequences'])
        ref = cast(str, compress_seq['reference'])
        vcf_metadata = cast(VCF_Metadata, compress_seq['metadata'])
    else:
        aln = correct_alignment(args.alignment, correct_nuc)
        ref = None
        vcf_metadata = None
        if args.root_sequence:
            for fmt in ['fasta', 'genbank']:
                try:
                    ref = str(read_single_sequence(args.root_sequence, format=fmt).seq).upper()
                    break
                except:
                    pass
            if ref is None:
                raise AugurError(f"could not read root sequence from {args.root_sequence}")
    if ref:
        ref = correct_nuc(ref)

    return (ref, aln, vcf_metadata)

def _to_ancestral_json(anc: Ancestral_Reconstruction) -> Ancestral_JSON:
    """Convert the results of a nucleotide ancestral reconstruction into the JSON
    format used by other augur tools"""
    root_seq = anc['root_seq']
    j: Ancestral_JSON = {
        'annotations': {
            'nuc': {'start': 1, 'end': len(root_seq), 'strand': '+', 'type': 'source'},
        },
        'reference': {'nuc': root_seq,},
        'nodes': anc['mutations']['nodes'],
    }
    if anc['mutations']['mask'] is not None:
        j['mask'] = "".join(['1' if x else '0' for x in anc["mutations"]["mask"]])
    return j

def _validate_translated_consistency(anc_seqs, aa_result, feat, gene, T) -> None:
    """Compare the independently reconstructed AA sequences with translations
    of the reconstructed nucleotide sequences for a given gene. Reports any
    inconsistencies as warnings.

    Parameters
    ----------
    anc_seqs : Ancestral_JSON
        Contains reconstructed nucleotide sequences for all nodes.
    aa_result : Ancestral_Reconstruction
        The independent AA reconstruction for this gene.
    feat : Bio.SeqFeature.SeqFeature
        The gene feature, used to extract the nucleotide region.
    gene : str
        Gene name, for reporting.
    T : Bio.Phylo.BaseTree.Tree
        The tree (used to iterate over all nodes).
    """
    nodes_with_differences = {'tips': 0, 'internal': 0}
    difference_counts: list[int] = []

    for node in T.find_clades():
        nuc_seq = anc_seqs['nodes'].get(node.name, {}).get('sequence')
        assert nuc_seq is not None
        nuc_translated = safe_translate(str(feat.extract(Seq(nuc_seq))))
        # re-reconstruct the inferred protein sequence since we don't store them on anc_seqs
        aa_seq = aa_result['tt'].sequence(node, as_string=True, reconstructed=True)
        assert len(aa_seq)==len(nuc_translated)
    
        if nuc_translated != aa_seq:
            nodes_with_differences["tips" if node.is_terminal() else "internal"] += 1
            difference_counts.append(sum(c1!=c2 for c1,c2 in zip(nuc_translated, aa_seq)))

    if len(difference_counts) > 0:
        difference_counts.sort()
        median = difference_counts[floor(len(difference_counts)/2)]
        print(
            f"WARNING: gene {gene!r} has differences between the independently reconstructed"
            f" amino acid sequence and a translation of the reconstructed nucleotide sequence."
            f" {nodes_with_differences['tips']} terminal nodes and"
            f" {nodes_with_differences['internal']} internal nodes were different."
            f" Median residue mismatch count: {median:,}"
            f" (range: {difference_counts[0]:,} - {difference_counts[-1]:,})",
            file=sys.stderr,
        )


def reconstruct_translations(
    anc_seqs: None|Ancestral_JSON,
    ref: str|None,
    T: Tree,
    genes_arg: list[str], # filename or list of genes
    annotation_fname: str,
    translations_fname_pattern: str,
    infer_ambiguous: bool,
    fill_overhangs: bool,
    marginal: bool,
    rng_seed: int,
    output_fname_pattern: str|None,
    report_inconsistent_translation: bool,
) -> Ancestral_JSON:
    genes = parse_genes_argument(genes_arg)
    if genes is None or not len(genes):
        raise AugurError("Empty list of genes provided")

    if len(genes)>1:
        if GENE_PATTERN not in translations_fname_pattern:
            raise AugurError(f"--translations must contain {GENE_PATTERN} for multiple-gene amino acid reconstructions")
        if output_fname_pattern and GENE_PATTERN not in output_fname_pattern:
            raise AugurError(f"--output-translations must contain {GENE_PATTERN} for multiple-gene amino acid reconstructions")

    ## load features; only requested features if genes given
    # Note: It's plausible to run a single-gene, no-nuc reconstruction without annotations if we want to allow that.
    from .io.sequences import load_features
    features = load_features(annotation_fname, genes)

    # anc_seqs was populated by the nucleotide reconstruction, which is optional. Create an ~empty structure if we
    # didn't run nucleotide reconstruction.
    if anc_seqs is None:
        anc_seqs = {
            "annotations": {},
            "reference": {},
            "nodes": {n.name: {} for n in T.find_clades()}
        }
    else:
        # Ensure the already-created nuc annotation coordinates match those parsed from the reference file
        assert 'nuc' in anc_seqs['annotations']
        if (features['nuc'].location.start+1 != anc_seqs['annotations']['nuc']['start'] or
            features['nuc'].location.end != anc_seqs['annotations']['nuc']['end']):
            raise AugurError(f"The 'nuc' annotation coordinates parsed from {annotation_fname!r} ({features['nuc'].location.start+1}..{features['nuc'].location.end})"
                f" don't match the provided sequence data coordinates ({anc_seqs['annotations']['nuc']['start']}..{anc_seqs['annotations']['nuc']['end']}).")

    correct_aa = _make_seq_corrector('aa')

    print("Read in {} features from reference sequence file".format(len(features)))
    for gene in genes:
        print(f"Processing gene: {gene}")
        fname = translations_fname_pattern.replace(GENE_PATTERN, gene) # GENE_PATTERN may not be in the filename if len(genes)==1
        feat = features[gene]
        reference_sequence = str(feat.extract(Seq(ref)).translate()) if ref else None
        if reference_sequence:
            reference_sequence = correct_aa(reference_sequence)

        aa_aln = correct_alignment(fname, correct_aa)
        aa_result = run_ancestral(T, aa_aln, reference_sequence=reference_sequence, is_vcf=False, fill_overhangs=fill_overhangs,
                                    marginal=marginal, infer_ambiguous=infer_ambiguous, alphabet='aa', rng_seed=rng_seed)
        len_translated_alignment = aa_result['tt'].data.full_length*3
        if len_translated_alignment != len(feat):
            raise AugurError(f"length of translated alignment for {gene} ({len_translated_alignment})"
                    f" does not match length of reference feature ({len(feat)})."
                    " Please make sure that the annotation matches the translations.")

        for key, node in anc_seqs['nodes'].items():
            if 'aa_muts' not in node:
                node['aa_muts'] = {}
            node['aa_muts'][gene] = aa_result['mutations']['nodes'][key]['muts']

            # Add amino acid sequences to the root node of the tree.
            if key == T.root.name:
                if "aa_sequences" not in node:
                    node["aa_sequences"] = {}

                node["aa_sequences"][gene] = aa_result['tt'].sequence(T.root, as_string=True, reconstructed=True)
        
        if report_inconsistent_translation:
            _validate_translated_consistency(anc_seqs, aa_result, feat, gene, T)

        anc_seqs['reference'][gene] = aa_result['root_seq']
        anc_seqs['annotations'].update(genome_features_to_auspice_annotation({gene: feat}, annotation_fname))

        # For each translated gene, save ancestral amino acid sequences to FASTA
        if output_fname_pattern:
            with open_file(output_fname_pattern.replace(GENE_PATTERN, gene), "w") as oh:
                for node in aa_result["tt"].tree.find_clades():
                    oh.write(f">{node.name}\n{aa_result['tt'].sequence(node, as_string=True, reconstructed=True)}\n")

    return anc_seqs


def run(args: argparse.Namespace):

    mode = validate_arguments(args)

    T = _read_tree(args.tree)

    import treetime
    print("\nInferred ancestral sequence states using TreeTime:"
          "\n\tSagulenko et al. TreeTime: Maximum-likelihood phylodynamic analysis"
          "\n\tVirus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731\n")
    print(f"augur ancestral is using TreeTime version {treetime.version}")

    # Infer ambiguous bases if the user has requested that we infer them (either
    # explicitly or by default) and the user has not explicitly requested that
    # we keep them.
    infer_ambiguous = bool(args.infer_ambiguous and not args.keep_ambiguous)
    rng_seed: int = args.seed
    fill_overhangs = not args.keep_overhangs
    marginal_inference = bool(args.inference=='marginal')
    
    # Nucleotide reconstruction is optional (we can reconstruct protein-only datasets)
    anc_seqs: None|Ancestral_JSON = None
    ref: str|None = None
    if mode.nuc_reconstruction:
        # read sequence data, correcting any invalid nucleotide states
        ref, aln, vcf_metadata = _read_sequence_data(args, mode.is_vcf)
        anc = run_ancestral(T, aln, reference_sequence=ref if ref else None, is_vcf=mode.is_vcf, fill_overhangs=fill_overhangs,
                            full_sequences=not mode.is_vcf, marginal=marginal_inference, infer_ambiguous=infer_ambiguous, alphabet='nuc',
                            rng_seed=rng_seed)
        anc_seqs = _to_ancestral_json(anc)

    if mode.aa_reconstruction:
        # For protein reconstruction by gene, read the already-translated (AA) FASTAs and reconstruct across the tree
        # Results are stored in the provided `anc_seqs` object and written to FASTA if requested
        assert not mode.is_vcf # guaranteed by validate_arguments() but good to double check
        anc_seqs = reconstruct_translations(anc_seqs, ref, T, args.genes, args.annotation, args.translations,
            infer_ambiguous, fill_overhangs, marginal_inference, rng_seed,
            args.output_translations, args.report_inconsistent_translation)
    
    default_json_fname = '.'.join(args.alignment.split('.')[:-1]) + '_mutations.json' if args.alignment else None
    out_name = get_json_name(args, default_json_fname)
    # use NodeDataObject to perform validation on the file before it's written
    NodeDataObject(anc_seqs, out_name, args.validation_mode)

    write_augur_json(anc_seqs, out_name)
    print("ancestral mutations written to", out_name, file=sys.stdout)

    if args.output_sequences:
        assert mode.nuc_reconstruction and not mode.is_vcf
        records = [
            SeqRecord(Seq(node_data["sequence"]), id=node_name, description="")
            for node_name, node_data in anc_seqs["nodes"].items()
        ]
        SeqIO.write(records, args.output_sequences, "fasta")
        print("ancestral sequences FASTA written to", args.output_sequences, file=sys.stdout)

    # output VCF including new ancestral seqs
    if args.output_vcf:
        assert mode.nuc_reconstruction and mode.is_vcf
        tree_dict = anc['tt'].get_tree_dict(keep_var_ambigs=not infer_ambiguous)
        tree_dict['metadata'] = vcf_metadata
        write_vcf(tree_dict, args.output_vcf, anc_seqs.get('mask'))
        print("Mutations, including for ancestral nodes, exported as VCF to", args.output_vcf, file=sys.stdout)

    return 0
