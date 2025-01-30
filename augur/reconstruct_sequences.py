"""
Reconstruct alignments from mutations inferred on the tree
"""

from Bio import SeqIO, Seq, SeqRecord, Phylo
from .io.file import open_file
from .io.sequences import read_sequences
from .utils import read_node_data



def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("reconstruct-sequences", help=__doc__)
    parser.add_argument('--tree', required=True, help="tree as Newick file")
    parser.add_argument('--gene', type=str, help="gene to translate (list or file containing list)")
    parser.add_argument('--mutations', required=True, type=str, help="json file containing mutations "
                            "mapped to each branch and the sequence of the root.")
    parser.add_argument('--vcf-aa-reference', type=str, help='fasta file of the reference gene translations for VCF format')
    parser.add_argument('--internal-nodes', action='store_true', help="include sequences of internal nodes in output")
    parser.add_argument('--output', type=str)
    return parser


def get_sequence(pseq, muts):
    """reconstruct a child sequence from that of its parent and the mutations

    Parameters
    ----------
    pseq : str
        sequence of the parent as string, won't be changed!
    muts : list
        mutations of the type K135T where K is the parental state, 135 is the
        position in 1-based numbering, and T is the new state

    Returns
    -------
    str
        reconstructed sequence
    """
    pseq_list = list(pseq)
    for mut in muts:
        new_state = mut[-1]
        pos = int(mut[1:-1])-1
        assert pseq_list[pos]==mut[0]
        pseq_list[pos]=new_state

    return "".join(pseq_list)


def load_alignments(sequence_files, gene_names):
    from Bio import AlignIO
    alignments = {}
    for fname, gene in zip(sequence_files, gene_names):
        alignments[gene] = AlignIO.read(fname, 'fasta')
    return alignments


def run(args):
    ## read tree and data, if reading data fails, return with error code
    tree = Phylo.read(args.tree, 'newick')

    #check whether VCF - by existance of VCF reference file
    is_vcf = args.vcf_aa_reference is not None

    ## check file format and read in sequences
    node_data = read_node_data(args.mutations, args.tree)
    if node_data is None:
        print("ERROR: could not read mutation data "+ ("(incl sequences)" if not is_vcf else ""))
        return 1

    root_node = tree.root.name

    #if VCF, read in the reference seq for each gene, put on root
    if(is_vcf):
        node_data["nodes"][root_node]['aa_sequences'] = {}
        with open_file(args.vcf_aa_reference) as handle:
            for record in read_sequences(handle, format="fasta"):
                if record.id==args.gene:
                    #'root' may not be same as 'reference', so apply any mutations at root here!
                    node_data["nodes"][root_node]['aa_sequences'][record.id] = get_sequence(str(record.seq), node_data["nodes"][root_node]["aa_muts"][record.id])

    # check that root node has sequences for each requested gene
    if "aa_sequences" not in node_data["nodes"][root_node]:
        print("ERROR: ancestral sequences are not provided")
        return 1
    if not args.gene in node_data["nodes"][root_node]['aa_sequences']:
        print("ERROR: ancestral sequences missing for gene",args.gene)
        return 1

    # gather all reconstructed sequences
    sequences = dict()
    is_terminal = {}

    sequences[root_node] = node_data["nodes"][root_node]['aa_sequences'][args.gene]
    is_terminal[root_node] = False
    for node in tree.get_nonterminals(order='preorder'):
        parent_sequence = sequences[node.name]
        for child in node:
            sequences[child.name] = get_sequence(parent_sequence, node_data["nodes"][child.name]["aa_muts"][args.gene])
            is_terminal[child.name] = child.is_terminal()

    # write alignments to file
    seqs = [SeqRecord.SeqRecord(seq=Seq.Seq(sequences[strain]), id=strain, name=strain, description='')
            for strain in sequences if is_terminal[strain] or args.internal_nodes]
    SeqIO.write(seqs, args.output, 'fasta')
