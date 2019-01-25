"""
Filter and subsample a sequence set.
"""

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import random, os, re
import numpy as np
from .utils import read_metadata, get_numerical_dates, run_shell_command

comment_char = '#'

def read_vcf(compressed, input_file):
    import gzip
    opn = gzip.open if compressed else open

    with opn(input_file, mode='rt') as f: #'rt' necessary for gzip
        for line in f:
            if line[0:2] == "#C":
                header = line.strip().split('\t')
                seq_keep = header[header.index("FORMAT")+1:]
                all_seq = seq_keep.copy() #because we need 'seqs to remove' for VCF
                return seq_keep, all_seq


def write_vcf(compressed, input_file, output_file, dropped_samps):
    #Read in/write out according to file ending
    inCall = "--gzvcf" if compressed else "--vcf"
    outCall = "| gzip -c" if output_file.lower().endswith('.gz') else ""

    toDrop = " ".join(["--remove-indv "+s for s in dropped_samps])
    call = ["vcftools", toDrop, inCall, input_file, "--recode --stdout", outCall, ">", output_file]

    print("Filtering samples using VCFTools with the call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    # remove vcftools log file
    try:
        os.remove('out.log')
    except OSError:
        pass

def read_priority_scores(fname):
    priorities = defaultdict(float)
    if not os.path.isfile(fname):
        print("ERROR: priority file %s doesn't exist"%fname)
        return priorities

    with open(fname) as pfile:
        for l in pfile:
            f = l.strip().split()
            try:
                priorities[f[0]] = float(f[1])
            except:
                print("ERROR: malformatted priority:",l)

    return priorities


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--metadata', required=True, help="metadata associated with sequences")
    parser.add_argument('--min-date', type=float, help="minimal cutoff for numerical date")
    parser.add_argument('--max-date', type=float, help="maximal cutoff for numerical date")
    parser.add_argument('--min-length', type=int, help="minimal length of the sequences")
    parser.add_argument('--non-nucleotide', action='store_true', help="exclude sequences that contain illegal characters")
    parser.add_argument('--exclude', type=str, help="file with list of strains that are to be excluded")
    parser.add_argument('--include', type=str, help="file with list of strains that are to be included regardless of priorities or subsampling")
    parser.add_argument('--priority', type=str, help="file with list priority scores for sequences (strain\tpriority)")
    parser.add_argument('--sequences-per-group', type=int, help="subsample to no more than this number of sequences per category")
    parser.add_argument('--group-by', nargs='+', help="categories with respect to subsample; two virtual fields, \"month\" and \"year\", are supported if they don't already exist as real fields but a \"date\" field does exist")
    parser.add_argument('--exclude-where', nargs='+',
                                help="Exclude samples matching these conditions. Ex: \"host=rat\" or \"host!=rat\". Multiple values are processed as OR (matching any of those specified will be excluded), not AND")
    parser.add_argument('--include-where', nargs='+',
                                help="Include samples with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. This rule is applied last and ensures any sequences matching these rules will be included.")
    parser.add_argument('--output', '-o', help="output file")


def run(args):
    '''
    filter and subsample a set of sequences into an analysis set
    '''
    #Set flags if VCF
    is_vcf = False
    is_compressed = False
    if any([args.sequences.lower().endswith(x) for x in ['.vcf', '.vcf.gz']]):
        is_vcf = True
        if args.sequences.lower().endswith('.gz'):
            is_compressed = True


    ####Read in files

    #If VCF, open and get sequence names
    if is_vcf:
        seq_keep, all_seq = read_vcf(is_compressed, args.sequences)

    #if Fasta, read in file to get sequence names and sequences
    else:
        seqs = {seq.id:seq for seq in SeqIO.parse(args.sequences, 'fasta')}
        seq_keep = list(seqs.keys())
        all_seq = seq_keep.copy()

    meta_dict, meta_columns = read_metadata(args.metadata)


    #####################################
    #Filtering steps
    #####################################

    # remove sequences without meta data
    tmp = [ ]
    for seq_name in seq_keep:
        if seq_name in meta_dict:
            tmp.append(seq_name)
        else:
            print("No meta data for %s, excluding from all further analysis."%seq_name)
    seq_keep = tmp

    # remove strains explicitly excluded by name
    # read list of strains to exclude from file and prune seq_keep
    if args.exclude and os.path.isfile(args.exclude):
        with open(args.exclude, 'r') as ifile:
            to_exclude = set()
            for line in ifile:
                if line[0] != comment_char:
                    # strip whitespace and remove all text following comment character
                    exclude_name = line.split(comment_char)[0].strip()
                    to_exclude.add(exclude_name)
        seq_keep = [seq_name for seq_name in seq_keep if seq_name not in to_exclude]

    # exclude strain my metadata field like 'host=camel'
    # match using lowercase
    if args.exclude_where:
        for ex in args.exclude_where:
            try:
                col, val = re.split(r'!?=', ex)
            except (ValueError,TypeError):
                print("invalid --exclude-where clause \"%s\", should be of from property=value or property!=value"%ex)
            else:
                to_exclude = set()
                for seq_name in seq_keep:
                    if "!=" in ex: # i.e. property!=value requested
                        if meta_dict[seq_name].get(col,'unknown').lower() != val.lower():
                            to_exclude.add(seq_name)
                    else: # i.e. property=value requested
                        if meta_dict[seq_name].get(col,'unknown').lower() == val.lower():
                            to_exclude.add(seq_name)
                seq_keep = [seq_name for seq_name in seq_keep if seq_name not in to_exclude]

    # filter by sequence length
    if args.min_length:
        if is_vcf: #doesn't make sense for VCF, ignore.
            print("WARNING: Cannot use min_length for VCF files. Ignoring...")
        else:
            seq_keep_by_length = []
            for seq_name in seq_keep:
                sequence = seqs[seq_name].seq
                length = sum(map(lambda x: sequence.count(x), ["a", "t", "g", "c", "A", "T", "G", "C"]))
                if length >= args.min_length:
                    seq_keep_by_length.append(seq_name)
            seq_keep = seq_keep_by_length

    # filter by date
    if (args.min_date or args.max_date) and 'date' in meta_columns:
        dates = get_numerical_dates(meta_dict, fmt="%Y-%m-%d")
        seq_keep = [s for s in seq_keep if dates[s] is not None]
        if args.min_date:
            seq_keep = [s for s in seq_keep if (np.isscalar(dates[s]) or all(dates[s])) and np.max(dates[s])>args.min_date]
        if args.max_date:
            seq_keep = [s for s in seq_keep if (np.isscalar(dates[s]) or all(dates[s])) and np.min(dates[s])<args.max_date]

    # exclude sequences with non-nucleotide characters
    if args.non_nucleotide:
        good_chars = {'A', 'C', 'G', 'T', '-', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'D', 'H', 'B', 'V', '?'}
        seq_keep = [s for s in seq_keep if len(set(str(seqs[s].seq).upper()).difference(good_chars))==0]

    # subsampling. This will sort sequences into groups by meta data fields
    # specified in --group-by and then take at most --sequences-per-group
    # from each group. Within each group, sequences are optionally sorted
    # by a priority score specified in a file --priority
    if args.group_by and args.sequences_per_group:
        spg = args.sequences_per_group
        seq_names_by_group = defaultdict(list)

        for seq_name in seq_keep:
            group = []
            m = meta_dict[seq_name]
            # collect group specifiers
            for c in args.group_by:
                if c in m:
                    group.append(m[c])
                elif c in ['month', 'year'] and 'date' in m:
                    try:
                        year = int(m["date"].split('-')[0])
                    except:
                        print("WARNING: no valid year, skipping",seq_name, m["date"])
                        continue
                    if c=='month':
                        try:
                            month = int(m["date"].split('-')[1])
                        except:
                            month = random.randint(1,12)
                        group.append((year, month))
                    else:
                        group.append(year)
                else:
                    group.append('unknown')
            seq_names_by_group[tuple(group)].append(seq_name)

        #If didnt find any categories specified, all seqs will be in 'unknown' - but don't sample this!
        if len(seq_names_by_group)==1 and 'unknown' in seq_names_by_group:
            print("WARNING: The specified group-by categories (%s) were not found."%args.group_by,
                  "No sequences-per-group sampling will be done.")
            if any([x in args.group_by for x in ['year','month']]):
                print("Note that using 'year' or 'year month' requires a column called 'date'.")
            print("\n")
        else:
            # Check to see if some categories are missing to warn the user
            group_by = set(['date' if cat in ['year','month'] else cat
                            for cat in args.group_by])
            missing_cats = [cat for cat in group_by if cat not in meta_columns]
            if missing_cats:
                print("WARNING:")
                if any([cat != 'date' for cat in missing_cats]):
                    print("\tSome of the specified group-by categories couldn't be found: ",
                          ", ".join([str(cat) for cat in missing_cats if cat != 'date']))
                if any([cat == 'date' for cat in missing_cats]):
                    print("\tA 'date' column could not be found to group-by year or month.")
                print("\tFiltering by group may behave differently than expected!\n")

            if args.priority: # read priorities
                priorities = read_priority_scores(args.priority)

            # subsample each groups, either by taking the spg highest priority strains or
            # sampling at random from the sequences in the group
            seq_subsample = []
            for group, sequences_in_group in seq_names_by_group.items():
                if args.priority: #sort descending by priority
                    seq_subsample.extend(sorted(sequences_in_group, key=lambda x:priorities[x], reverse=True)[:spg])
                else:
                    seq_subsample.extend(sequences_in_group if len(sequences_in_group)<=spg
                                         else random.sample(sequences_in_group, spg))

            seq_keep = seq_subsample

    # force include sequences specified in file.
    # Note that this might re-add previously excluded sequences
    # Note that we are also not checking for existing meta data here
    if args.include and os.path.isfile(args.include):
        with open(args.include, 'r') as ifile:
            to_include = set([line.strip() for line in ifile if line[0]!=comment_char])

        for s in to_include:
            if s not in seq_keep:
                seq_keep.append(s)

    # add sequences with particular meta data attributes
    if args.include_where:
        to_include = []
        for ex in args.include_where:
            try:
                col, val = ex.split("=")
            except (ValueError,TypeError):
                print("invalid include clause %s, should be of from property=value"%ex)
                continue

            # loop over all sequences and re-add sequences
            for seq_name in all_seq:
                if seq_name in meta_dict:
                    if meta_dict[seq_name].get(col)==val:
                        to_include.append(seq_name)
                else:
                    print("WARNING: no metadata for %s, skipping"%seq_name)
                    continue

        for s in to_include:
            if s not in seq_keep:
                seq_keep.append(s)

    ####Write out files

    if is_vcf:
        #get the samples to be deleted, not to keep, for VCF
        dropped_samps = list(set(all_seq) - set(seq_keep))
        if len(dropped_samps) == len(all_seq): #All samples have been dropped! Stop run, warn user.
            print("ERROR: All samples have been dropped! Check filter rules and metadata file format.")
            return 1
        write_vcf(is_compressed, args.sequences, args.output, dropped_samps)

    else:
        seq_to_keep = [seq for id,seq in seqs.items() if id in seq_keep]
        if len(seq_to_keep) == 0:
            print("ERROR: All samples have been dropped! Check filter rules and metadata file format.")
            return 1
        SeqIO.write(seq_to_keep, args.output, 'fasta')
