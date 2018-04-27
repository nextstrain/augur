import sys,argparse

def run(args):
    parser = argparse.ArgumentParser("Translate the nucleotide alignments")
    parser.add_argument('--tree', required=True, help="tree file")
    parser.add_argument('--metadata', required=True, help="file with meta data")
    parser.add_argument('--attribute', default='region',
                        help='meta data field to perform discrete reconstruction on')
    parser.add_argument('--confidence',action="store_true",
                        help='record the distribution of subleading mugration states')
    parser.add_argument('--vcf', action='store_true', default=False,
                        help="sequence is in VCF format")
    args = parser.parse_args(args)


if __name__ == '__main__':
    run(sys.argv)