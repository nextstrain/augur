"""
Translate gene regions from nucleotides to amino acids.

Translates nucleotide sequences of nodes in a tree to amino acids for gene
regions of the annotated features of the provided reference sequence.
Each node then gets assigned a list of amino acid mutations for any position
that has a mismatch between its own amino acid sequence and its parent's sequence.
The reference amino acid sequences, genome annotations, and node amino acid
mutations are output to a node-data JSON file.

.. note::

    The mutation positions in the node-data JSON are one-based.
"""

import sys
import numpy as np
from Bio import SeqIO, Seq, SeqRecord, Phylo
from .io.vcf import write_VCF_translation, is_vcf as is_filename_vcf
from .utils import parse_genes_argument, read_node_data, load_features, \
    write_json, get_json_name, genome_features_to_auspice_annotation
from treetime.vcf_utils import read_vcf
from augur.errors import AugurError
from textwrap import dedent
from .types import ValidationMode
from .util_support.node_data_file import NodeDataObject
from .export_v2 import validation_mode_help_message

class MissingNodeError(Exception):
    pass

class MismatchNodeError(Exception):
    pass

class NoVariationError(Exception):
    pass

def safe_translate(sequence, report_exceptions=False):
    """Returns an amino acid translation of the given nucleotide sequence accounting
    for gaps in the given sequence.

    Optionally, returns a tuple of the translated sequence and whether an
    exception was raised during initial translation.

    Examples
    --------
    >>> safe_translate("ATG")
    'M'
    >>> safe_translate("ATGGT-")
    'MX'
    >>> safe_translate("ATG---")
    'M-'
    >>> safe_translate("ATGTAG")
    'M*'
    >>> safe_translate("")
    ''
    >>> safe_translate("ATGT")
    'MX'
    >>> safe_translate("ATG", report_exceptions=True)
    ('M', False)
    >>> safe_translate("ATGA-G", report_exceptions=True)
    ('MX', True)
    """
    from Bio.Data.CodonTable import TranslationError
    from Bio.Seq import CodonTable
    translation_exception = False

    #sequences not mod 3 give messy BiopythonWarning, so avoid by padding.
    if len(sequence)%3:
        sequence_padded = sequence + "N"*(3-len(sequence)%3)
    else:
        sequence_padded = sequence
    try:
        # Attempt translation by extracting the sequence according to the
        # BioPhython SeqFeature in frame gaps of three will translate as '-'
        translated_sequence = str(Seq.Seq(sequence_padded).translate(gap='-'))
    except TranslationError:
        translation_exception = True
        # Any other codon like '-AA' or 'NNT' etc will fail. Translate codons
        # one by one.
        codon_table  = CodonTable.ambiguous_dna_by_name['Standard'].forward_table
        str_seq = str(sequence_padded)
        codons = np.frombuffer(str_seq[:len(str_seq) - len(str_seq) % 3].encode(), dtype='S3').astype("U")
        assert len(codons) > 0
        aas = []

        for c in codons:
            # Parse result of single codon translation, add amino acids as
            # appropriate.
            try:
                aa = codon_table.get(c)
                if aa is None:
                    if c == '---':
                        aas.append('-')
                    else:
                        aas.append('X')
                else:
                    aas.append(aa)
            except (TranslationError, ValueError):
                aas.append('X')

        translated_sequence = "".join(aas)

    if report_exceptions:
        return translated_sequence, translation_exception
    else:
        return translated_sequence


def translate_feature(aln, feature):
    '''
    Translates a subsequence of input nucleotide sequences.

    Parameters
    ----------
    aln : dict
        sequences indexed by node name

    feature : Bio.Seq.Seq
        BioPython sequence feature

    Returns
    -------
    dict :
        translated sequences indexed by node name

    '''
    translations = {}
    for sname, seq in aln.items():
        aa_seq = safe_translate(str(feature.extract(seq)))
        translations[sname] = aa_seq

    return translations


def translate_vcf_feature(sequences, ref, feature, feature_name):
    '''Translates a subsequence of input nucleotide sequences.

    Parameters
    ----------
    sequences : dict
        TreeTime format dictionary from VCF-input of sequences indexed by node name

    ref :
        reference alignment the VCF was mapped to

    feature : Bio.Seq.Seq
        BioPython sequence feature

    Returns
    -------
    dict :
        translated reference gene, positions of AA differences, and AA
        differences indexed by node name

    Raises
    ------
    NoVariationError : if no variable sites within this feature (across all sequences)
    '''
    def str_reverse_comp(str_seq):
        #gets reverse-compliment of a string and returns it as a string
        seq_str = Seq(str_seq)
        return str(seq_str.reverse_complement())

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    prot = {}
    prot['sequences'] = {}
    prot['positions'] = []

    refNuc = str(feature.extract( SeqRecord(seq=Seq(ref)) ).seq)
    # Need to get ref translation to store. check if multiple of 3 for sanity.
    # will be padded in safe_translate if not
    if len(refNuc)%3:
        print(f"Gene length of {feature_name!r} is not a multiple of 3. will pad with N", file=sys.stderr)

    ref_aa_seq = safe_translate(refNuc)
    prot['reference'] = ref_aa_seq

    start = int(feature.location.start) #1 less than GFF file; Bio 'pythons' the numbering
    end = int(feature.location.end)

    for seqk in sequences.keys():
        varSite = np.array(list(sequences[seqk].keys()))  #get positions where nuc diffs
        #reduce to only those within current gene
        geneVarSites = np.logical_and(varSite >= start, varSite < end)
        nucVarSites = varSite[geneVarSites] #translate this back to nuc position
        genNucSites = nucVarSites-start #get it in position within the gene to ensure in frame

        #Translate just the codon this nuc diff is in, and find out which AA loc
        #But need numbering to be w/in protin, not whole genome
        if feature.location.strand == -1:
            aaRepLocs = {(end-start-i-1)//3:safe_translate( str_reverse_comp( "".join([sequences[seqk][key+start]
                                    if key+start in sequences[seqk].keys() else ref[key+start]
                                for key in range(i-i%3,i+3-i%3)]) ))
                            for i in genNucSites}
        else:
            aaRepLocs = {i//3:safe_translate( "".join([sequences[seqk][key+start]
                                if key+start in sequences[seqk].keys() else ref[key+start]
                            for key in range(i-i%3,i+3-i%3)]) )
                        for i in genNucSites}

        aaRepLocsFinal = {}
        #remove if is a synonymous mutation!
        for key,val in aaRepLocs.items():
            if ref_aa_seq[key] != val:
                aaRepLocsFinal[key] = val
        aaRepLocs = aaRepLocsFinal

        #store the dict of differences
        prot['sequences'][seqk] = aaRepLocs
        #add to list of positions if needed
        for key in aaRepLocs.keys():
            if key not in prot['positions']:
                prot['positions'].append(key)

    prot['positions'].sort()

    # raise an error if no variable sites observed
    if len(prot['positions']) == 0:
        raise NoVariationError()
    return prot

def construct_mut(start, pos, end):
    return str(start) + str(pos) + str(end)

def assign_aa_vcf(tree, translations):
    aa_muts = {}

    #get mutations on the root
    root = tree.root
    aa_muts[root.name]={"aa_muts":{}}
    #If has no root node name, exit with error
    if root.name is None:
        print("\n*** Can't find node name for the tree root!")
        raise MissingNodeError()

    for fname, prot in translations.items():
        if root.name not in prot['sequences']:
            print("\n*** Can't find %s in the alignment provided!"%(root.name))
            raise MismatchNodeError()
        root_muts = prot['sequences'][root.name]
        tmp = []
        for pos in prot['positions']:
            if pos in root_muts:
                tmp.append(construct_mut(prot['reference'][pos], int(pos+1), root_muts[pos]))
        aa_muts[root.name]["aa_muts"][fname] = tmp

    for n in tree.get_nonterminals():
        for c in n:
            aa_muts[c.name]={"aa_muts":{}}
        for fname, prot in translations.items():
            if n.name not in prot['sequences']:
                print("\n*** Can't find %s in the alignment provided!"%(root.name))
                raise MismatchNodeError()
            n_muts = prot['sequences'][n.name]
            for c in n:
                tmp = []
                if c.name is None:
                    print("\n*** Internal node missing a node name!")
                    raise MissingNodeError()
                c_muts = prot['sequences'][c.name]
                for pos in prot['positions']:
                    #if pos in both, check if same
                    if pos in n_muts and pos in c_muts:
                        if n_muts[pos] != c_muts[pos]:
                            tmp.append(construct_mut(n_muts[pos], int(pos+1), c_muts[pos]))
                    elif pos in n_muts:
                        tmp.append(construct_mut(n_muts[pos], int(pos+1), prot['reference'][pos]))
                    elif pos in c_muts:
                        tmp.append(construct_mut(prot['reference'][pos], int(pos+1), c_muts[pos]))

                aa_muts[c.name]["aa_muts"][fname] = tmp

    aa_muts[root.name]['aa_sequences'] = {}
    for gene_name, gene_data in translations.items():
        root_seq = list(gene_data['reference'])
        for pos,alt in gene_data['sequences'][root.name].items():
            # pos is 0-based, <class 'numpy.int64'>
            root_seq[pos] = alt
        aa_muts[root.name]['aa_sequences'][gene_name] = "".join(root_seq)

    return aa_muts

def assign_aa_fasta(tree, translations, reference_translations):
    aa_muts = {}

    # Depending on how `augur ancestral` was run, the input JSON (nt_muts.json)
    # may or may not have mutations defined on the root node. Note that the 'muts'
    # array will always be present, but it can only contain mutations if a
    # --root-sequence was provided to augur ancestral.

    root = tree.get_nonterminals()[0]

    for n in tree.get_nonterminals():
        if n.name is None:
            print("\n*** Tree is missing node names!")
            raise MissingNodeError()
        for c in n:
            aa_muts[c.name]={"aa_muts":{}}
        for fname, aln in translations.items():
            for c in n:
                if c.name in aln and n.name in aln:
                    tmp = [construct_mut(a, int(pos+1), d) for pos, (a,d) in
                            enumerate(zip(aln[n.name], aln[c.name])) if a!=d]
                    aa_muts[c.name]["aa_muts"][fname] = tmp
                elif c.name not in aln and n.name not in aln:
                    print("\n*** Can't find %s OR %s in the alignment provided!"%(c.name, n.name))
                    raise MismatchNodeError()
                else:
                    print("no sequence pair for nodes %s-%s"%(c.name, n.name))

        if n==tree.root:
            # The aa_sequences at the root are simply the translation from the root-node input data
            aa_muts[n.name]={"aa_sequences":{}}
            for fname, aln in translations.items():
                if n.name in aln:
                    aa_muts[n.name]["aa_sequences"][fname] = "".join(aln[n.name])
            # The aa_muts are found by comparing the aa_sequence with the reference sequence
            aa_muts[n.name]['aa_muts'] = {}
            for fname, aln in translations.items():
                ref = reference_translations[fname]
                muts = [construct_mut(a, int(pos+1), d) for pos, (a,d) in enumerate(zip(ref, aln[n.name])) if a!=d]
                aa_muts[n.name]["aa_muts"][fname] = muts

    return aa_muts

def sequences_vcf(reference_fasta, vcf):
    """
    Extract the nucleotide variation in the VCF
    Returns a tuple
    [0] The sequences as a dict of dicts. sequences → <NODE_NAME> → <POS> → <ALT_NUC> where <POS> is a 0-based int
    [1] The sequence of the provided `reference_fasta` (string)
    """
    assert reference_fasta is not None
    compress_seq = read_vcf(vcf, reference_fasta)
    sequences = compress_seq['sequences']
    ref = compress_seq['reference']
    return (sequences, ref)

def sequences_json(node_data_json, tree, validation_mode):
    """
    Extract the full nuc sequence for each node in the provided node-data JSON.
    Returns a dict, keys are node names and values are a string of the genome sequence (nuc)
    """
    node_data = read_node_data(node_data_json, validation_mode=validation_mode)
    if node_data is None:
        raise AugurError("could not read node data (incl sequences)")
    # extract sequences from node meta data
    sequences = {}
    for k,v in node_data['nodes'].items():
        if 'sequence' in v:
            sequences[k] = v['sequence']
    tree_nodes = {c.name for c in tree.find_clades()}
    tree_nodes_missing_sequences = tree_nodes - set(sequences.keys())
    if len(tree_nodes_missing_sequences):
        raise AugurError(dedent(f"""\
            {len(tree_nodes_missing_sequences)} nodes on the tree are missing nucleotide sequences in the node-data JSON.
            These must be present under 'nodes' → <node_name> → 'sequence'.
            This error may originate from using 'augur ancestral' with VCF input; in this case try using VCF output from that command here.
            """))
    reference = node_data['reference']['nuc']
    return (reference, sequences)

def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("translate", help=__doc__)
    parser.add_argument('--tree', required=True, help="prebuilt Newick -- no tree will be built if provided")
    parser.add_argument('--ancestral-sequences', required=True, type=str, help='JSON (fasta input) or VCF (VCF input) containing ancestral and tip sequences')
    parser.add_argument('--reference-sequence', required=True,
                        help='GenBank or GFF file containing the annotation')
    parser.add_argument('--genes', nargs='+', action='extend', help="genes to translate (list or file containing list)")
    parser.add_argument('--output-node-data', type=str, help='name of JSON file to save aa-mutations to')
    parser.add_argument('--alignment-output', type=str, help="write out translated gene alignments. "
                                   "If a VCF-input, a .vcf or .vcf.gz will be output here (depending on file ending). If fasta-input, specify the file name "
                                   "like so: 'my_alignment_%%GENE.fasta', where '%%GENE' will be replaced by the name of the gene")
    parser.add_argument('--validation-mode', type=ValidationMode, choices=[mode for mode in ValidationMode], default=ValidationMode.ERROR, help=validation_mode_help_message)

    vcf_only = parser.add_argument_group(
        title="VCF specific",
        description="These arguments are only applicable if the input (--ancestral-sequences) is in VCF format."
    )
    vcf_only.add_argument('--vcf-reference', type=str, help='fasta file of the sequence the VCF was mapped to')
    vcf_only.add_argument('--vcf-reference-output', type=str, help="fasta file where reference sequence translations for VCF input will be written")

    return parser

def check_arg_combinations(args, is_vcf):
    """
    Check that provided arguments are compatible.
    Where possible we use argparse built-ins, but they don't cover everything we want to check.
    This checking shouldn't be used by downstream code to assume arguments exist, however by checking for
    invalid combinations up-front we can exit quickly.
    """

    if is_vcf:
        if not args.vcf_reference:
            raise AugurError("A reference FASTA (--vcf-reference) is required with VCF-format input")
    else:
        if args.vcf_reference or args.vcf_reference_output:
            raise AugurError("Arguments '--vcf-reference' and/or '--vcf-reference-output' are only applicable if the input ('--ancestral-sequences') is VCF")
    
    if args.alignment_output:
        if is_vcf:
            if not is_filename_vcf(args.alignment_output):
                raise AugurError("When using a VCF input the --alignment-output filename must also be a VCF file")
            if not args.vcf_reference_output:
                raise AugurError("When using a VCF input and --alignment-output, we now require you to specify the --vcf-reference-output as well")
        else:
            if is_filename_vcf(args.alignment_output):
                raise AugurError("When using a non-VCF input the --alignment-output filename must not be a VCF file")
    if args.vcf_reference_output and not args.alignment_output:
        raise AugurError("The VCF reference output (--vcf-reference-output) needs --alignment-output")


def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')
    is_vcf = is_filename_vcf(args.ancestral_sequences)
    check_arg_combinations(args, is_vcf)

    genes = parse_genes_argument(args.genes)

    ## load features; only requested features if genes given
    features = load_features(args.reference_sequence, genes)
    print("Read in {} features from reference sequence file".format(len(features)))

    ## Read in sequences & for each sequence translate each feature _except for_ the 'nuc' feature name
    ## Note that except for the 'nuc' annotation, `load_features` _only_ looks for 'gene' (GFF files) or 'CDS' (GenBank files)
    ## The reference translations are straight translations of the "reference.nuc" sequence in the JSON
    ## (see <https://github.com/nextstrain/augur/issues/1362> for more details here), _or_ for VCF input a translation
    ## from the provided FASTA reference
    translations = {}
    reference_translations = {}
    if is_vcf:
        (sequences, ref) = sequences_vcf(args.vcf_reference, args.ancestral_sequences)
        features_without_variation = []
        for fname, feat in features.items():
            if fname=='nuc':
                continue
            try:
                translations[fname] = translate_vcf_feature(sequences, ref, feat, fname)
                reference_translations[fname] = translations[fname]['reference']
            except NoVariationError:
                features_without_variation.append(fname)
        if len(features_without_variation):
            print("{} genes had no mutations and so have been be excluded.".format(len(features_without_variation)))  
    else:
        (reference, sequences) = sequences_json(args.ancestral_sequences, tree, args.validation_mode)
        translations = {fname: translate_feature(sequences, feat) for fname, feat in features.items() if fname!='nuc'}
        for fname, feat in features.items():
            if fname=='nuc':
                continue
            reference_translations[fname] = safe_translate(str(feat.extract(reference)))

    annotations = genome_features_to_auspice_annotation(features, args.reference_sequence, assert_nuc=True)

    ## determine amino acid mutations for each node
    try:
        if is_vcf:
            aa_muts = assign_aa_vcf(tree, translations)
        else:
            aa_muts = assign_aa_fasta(tree, translations, reference_translations)
    except MissingNodeError as err:
        print("\n*** ERROR: Some/all nodes have no node names!")
        print("*** Please check you are providing the tree output by 'augur refine'.")
        print("*** If you haven't run 'augur refine', please add node names to your tree by running:")
        print("*** augur refine --tree %s --output-tree <filename>.nwk"%(args.tree) )
        print("*** And use <filename>.nwk as the tree when running 'ancestral', 'translate', and 'traits'")
        return 1
    except MismatchNodeError as err:
        print("\n*** ERROR: Mismatch between node names in %s and in %s"%(args.tree, args.ancestral_sequences))
        print("*** Ensure you are using the same tree you used to run 'ancestral' as input here.")
        print("*** Or, re-run 'ancestral' using %s, then use the new %s as input here."%(args.tree, args.ancestral_sequences))
        return 1

    output_data = {'annotations':annotations, 'nodes':aa_muts, 'reference': reference_translations}
    out_name = get_json_name(args, '.'.join(args.tree.split('.')[:-1]) + '_aa-mutations.json')
    # use NodeDataObject to perform validation on the file before it's written
    NodeDataObject(output_data, out_name, args.validation_mode)

    write_json(output_data, out_name)
    print("amino acid mutations written to", out_name, file=sys.stdout)

    ## write alignments to file if requested
    if args.alignment_output:
        if is_vcf:
            assert is_filename_vcf(args.alignment_output)
            assert args.vcf_reference_output is not None
            write_VCF_translation(translations, args.alignment_output, args.vcf_reference_output)
        else:
            ## write fasta-style output if requested
            if '%GENE' in args.alignment_output:
                for fname, seqs in translations.items():
                    SeqIO.write([SeqRecord.SeqRecord(seq=Seq.Seq(s), id=sname, name=sname, description='')
                                 for sname, s in seqs.items()],
                                 args.alignment_output.replace('%GENE', fname), 'fasta')
            else:
                print("ERROR: alignment output file does not contain '%GENE', so will not be written.")
