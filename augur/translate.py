import sys,argparse

def run(args):
    parser = argparse.ArgumentParser("Translate the nucleotide alignments")
    parser.add_argument('-s', required=True, help="nucleotide alignment")
    parser.add_argument('--reference', required=True,
                        help='genbank file containing the annotation')
    parser.add_argument('--genes', nargs='+', help="genes to translate")
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")
    parser.add_argument('--assignMuts', action='store_true', default=False,
                        help="write amino acid mutations onto the tree")
    args = parser.parse_args()

    path = args.path


if __name__ == '__main__':
    run(sys.argv)