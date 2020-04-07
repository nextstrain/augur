"""
Filter and subsample a sequence set.
"""

from Bio import SeqIO
import pandas as pd
from collections import defaultdict
import random, os, re
import numpy as np
import sys
from .utils import read_metadata, get_numerical_dates, run_shell_command, shquote

comment_char = '#'


def read_vcf(filename):
    if filename.lower().endswith(".gz"):
        import gzip
        file = gzip.open(filename, mode="rt")
    else:
        file = open(filename)

    chrom_line = next(line for line in file if line.startswith("#C"))
    file.close()
    headers = chrom_line.strip().split("\t")
    sequences = headers[headers.index("FORMAT") + 1:]

    # because we need 'seqs to remove' for VCF
    return sequences, sequences.copy()


def write_vcf(input_filename, output_filename, dropped_samps):
    if _filename_gz(input_filename):
        input_arg = "--gzvcf"
    else:
        input_arg = "--vcf"

    if _filename_gz(output_filename):
        output_pipe = "| gzip -c"
    else:
        output_pipe = ""

    drop_args = ["--remove-indv " + shquote(s) for s in dropped_samps]

    call = ["vcftools"] + drop_args + [input_arg, shquote(input_filename), "--recode --stdout", output_pipe, ">", shquote(output_filename)]

    print("Filtering samples using VCFTools with the call:")
    print(" ".join(call))
    run_shell_command(" ".join(call), raise_errors = True)
    # remove vcftools log file
    try:
        os.remove('out.log')
    except OSError:
        pass

def read_priority_scores(fname):
    try:
        with open(fname) as pfile:
            return defaultdict(float, {
                elems[0]: float(elems[1])
                for elems in (line.strip().split() for line in pfile.readlines())
            })
    except Exception as e:
        print(f"ERROR: missing or malformed priority scores file {fname}", file=sys.stderr)
        raise e


def register_arguments(parser):
    parser.add_argument('--sequences', '-s', required=True, help="sequences in fasta or VCF format")
    parser.add_argument('--metadata', required=True, help="metadata associated with sequences")
    parser.add_argument('--min-date', type=float, help="minimal cutoff for numerical date")
    parser.add_argument('--max-date', type=float, help="maximal cutoff for numerical date")
    parser.add_argument('--min-length', type=int, help="minimal length of the sequences")
    parser.add_argument('--non-nucleotide', action='store_true', help="exclude sequences that contain illegal characters")
    parser.add_argument('--exclude', type=str, help="file with list of strains that are to be excluded")
    parser.add_argument('--include', type=str, help="file with list of strains that are to be included regardless of priorities or subsampling")
    parser.add_argument('--priority', type=str, help="file with list of priority scores for sequences (strain\tpriority)")
    parser.add_argument('--sequences-per-group', type=int, help="subsample to no more than this number of sequences per category")
    parser.add_argument('--group-by', nargs='+', help="categories with respect to subsample; two virtual fields, \"month\" and \"year\", are supported if they don't already exist as real fields but a \"date\" field does exist")
    parser.add_argument('--subsample-seed', help="random number generator seed to allow reproducible sub-sampling (with same input data). Can be number or string.")
    parser.add_argument('--exclude-where', nargs='+',
                                help="Exclude samples matching these conditions. Ex: \"host=rat\" or \"host!=rat\". Multiple values are processed as OR (matching any of those specified will be excluded), not AND")
    parser.add_argument('--include-where', nargs='+',
                                help="Include samples with these values. ex: host=rat. Multiple values are processed as OR (having any of those specified will be included), not AND. This rule is applied last and ensures any sequences matching these rules will be included.")
    parser.add_argument('--output', '-o', help="output file", required=True)


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

    ### Check users has vcftools. If they don't, a one-blank-line file is created which
    #   allows next step to run but error very badly.
    if is_vcf:
        from shutil import which
        if which("vcftools") is None:
            print("ERROR: 'vcftools' is not installed! This is required for VCF data. "
                  "Please see the augur install instructions to install it.")
            return 1

    ####Read in files

    #If VCF, open and get sequence names
    if is_vcf:
        seq_keep, all_seq = read_vcf(args.sequences)

    #if Fasta, read in file to get sequence names and sequences
    else:
        try:
            seqs = SeqIO.to_dict(SeqIO.parse(args.sequences, 'fasta'))
        except ValueError as error:
            print("ERROR: Problem reading in {}:".format(args.sequences))
            print(error)
            return 1
        seq_keep = list(seqs.keys())
        all_seq = seq_keep.copy()

    try:
        meta_dict, meta_columns = read_metadata(args.metadata)
    except ValueError as error:
        print("ERROR: Problem reading in {}:".format(args.metadata))
        print(error)
        return 1


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
    num_excluded_by_name = 0
    if args.exclude:
        try:
            with open(args.exclude, 'r') as ifile:
                to_exclude = set()
                for line in ifile:
                    if line[0] != comment_char:
                        # strip whitespace and remove all text following comment character
                        exclude_name = line.split(comment_char)[0].strip()
                        to_exclude.add(exclude_name)
            tmp = [seq_name for seq_name in seq_keep if seq_name not in to_exclude]
            num_excluded_by_name = len(seq_keep) - len(tmp)
            seq_keep = tmp
        except FileNotFoundError as e:
            print("ERROR: Could not open file of excluded strains '%s'" % args.exclude, file=sys.stderr)
            sys.exit(1)

    # exclude strain my metadata field like 'host=camel'
    # match using lowercase
    num_excluded_by_metadata = {}
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
                tmp = [seq_name for seq_name in seq_keep if seq_name not in to_exclude]
                num_excluded_by_metadata[ex] = len(seq_keep) - len(tmp)
                seq_keep = tmp

    # filter by sequence length
    num_excluded_by_length = 0
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
            num_excluded_by_length = len(seq_keep) - len(seq_keep_by_length)
            seq_keep = seq_keep_by_length

    # filter by date
    num_excluded_by_date = 0
    if (args.min_date or args.max_date) and 'date' in meta_columns:
        dates = get_numerical_dates(meta_dict, fmt="%Y-%m-%d")
        tmp = [s for s in seq_keep if dates[s] is not None]
        if args.min_date:
            tmp = [s for s in tmp if (np.isscalar(dates[s]) or all(dates[s])) and np.max(dates[s])>args.min_date]
        if args.max_date:
            tmp = [s for s in tmp if (np.isscalar(dates[s]) or all(dates[s])) and np.min(dates[s])<args.max_date]
        num_excluded_by_date = len(seq_keep) - len(tmp)
        seq_keep = tmp

    # exclude sequences with non-nucleotide characters
    num_excluded_by_nuc = 0
    if args.non_nucleotide:
        good_chars = {'A', 'C', 'G', 'T', '-', 'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'D', 'H', 'B', 'V', '?'}
        tmp = [s for s in seq_keep if len(set(str(seqs[s].seq).upper()).difference(good_chars))==0]
        num_excluded_by_nuc = len(seq_keep) - len(tmp)
        seq_keep = tmp

    # subsampling. This will sort sequences into groups by meta data fields
    # specified in --group-by and then take at most --sequences-per-group
    # from each group. Within each group, sequences are optionally sorted
    # by a priority score specified in a file --priority
    # Fix seed for the RNG if specified
    if args.subsample_seed:
        random.seed(args.subsample_seed)
    num_excluded_subsamp = 0
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
        if len(seq_names_by_group)==1 and ('unknown' in seq_names_by_group or ('unknown',) in seq_names_by_group):
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

            num_excluded_subsamp = len(seq_keep) - len(seq_subsample)
            seq_keep = seq_subsample

    # force include sequences specified in file.
    # Note that this might re-add previously excluded sequences
    # Note that we are also not checking for existing meta data here
    num_included_by_name = 0
    if args.include and os.path.isfile(args.include):
        with open(args.include, 'r') as ifile:
            to_include = set(
                [
                    line.strip()
                    for line in ifile
                    if line[0]!=comment_char and len(line.strip()) > 0
                ]
            )

        for s in to_include:
            if s not in seq_keep:
                seq_keep.append(s)
                num_included_by_name += 1

    # add sequences with particular meta data attributes
    num_included_by_metadata = 0
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
                num_included_by_metadata += 1

    ####Write out files

    if is_vcf:
        #get the samples to be deleted, not to keep, for VCF
        dropped_samps = list(set(all_seq) - set(seq_keep))
        if len(dropped_samps) == len(all_seq): #All samples have been dropped! Stop run, warn user.
            print("ERROR: All samples have been dropped! Check filter rules and metadata file format.")
            return 1
        write_vcf(args.sequences, args.output, dropped_samps)

    else:
        seq_to_keep = [seq for id,seq in seqs.items() if id in seq_keep]
        if len(seq_to_keep) == 0:
            print("ERROR: All samples have been dropped! Check filter rules and metadata file format.")
            return 1
        SeqIO.write(seq_to_keep, args.output, 'fasta')

    print("\n%i sequences were dropped during filtering" % (len(all_seq) - len(seq_keep),))
    if args.exclude:
        print("\t%i of these were dropped because they were in %s" % (num_excluded_by_name, args.exclude))
    if args.exclude_where:
        for key,val in num_excluded_by_metadata.items():
            print("\t%i of these were dropped because of '%s'" % (val, key))
    if args.min_length:
        print("\t%i of these were dropped because they were shorter than minimum length of %sbp" % (num_excluded_by_length, args.min_length))
    if (args.min_date or args.max_date) and 'date' in meta_columns:
        print("\t%i of these were dropped because of their date (or lack of date)" % (num_excluded_by_date))
    if args.non_nucleotide:
        print("\t%i of these were dropped because they had non-nucleotide characters" % (num_excluded_by_nuc))
    if args.group_by and args.sequences_per_group:
        seed_txt = ", using seed {}".format(args.subsample_seed) if args.subsample_seed else ""
        print("\t%i of these were dropped because of subsampling criteria%s" % (num_excluded_subsamp, seed_txt))

    if args.include and os.path.isfile(args.include):
        print("\n\t%i sequences were added back because they were in %s" % (num_included_by_name, args.include))
    if args.include_where:
        print("\t%i sequences were added back because of '%s'" % (num_included_by_metadata, args.include_where))

    print("%i sequences have been written out to %s" % (len(seq_keep), args.output))


def _filename_gz(filename):
    return filename.lower().endswith(".gz")
