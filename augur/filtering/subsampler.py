from collections import defaultdict
import functools
import random


class Subsampler:
    def __init__(self, *, priorities_fname, group_by, sequences_per_group, seed=None):
        self.priorities_fname = priorities_fname
        self.discriminators = group_by
        self.sequences_per_group = sequences_per_group
        self.seed = seed

    def is_affected(self, sequences):
        """
        Subsampling sorts sequences into groups by metadata fields specified in --group-by and then take at
        most --sequences-per-group from each group. Within each group, sequences are optionally sorted by a
        priority score specified in a file --priority.
        The random seed may be specified by --subsample-seed.
        """
        return False  # TODO
        if self.seed:
            random.seed(self.seed)

        seq_names_by_group = self.divide_into_groups(sequences)

        # If didnt find any categories specified, all seqs will be in 'unknown' - but don't sample this!
        if len(seq_names_by_group) == 1 and ("unknown" in seq_names_by_group or ("unknown",) in seq_names_by_group):
            print(f"WARNING: The specified group-by categories ({self.discriminators}) were not found. No sequences-per-group sampling will be done.")
            if any([x in self.discriminators for x in ['year', 'month']]):
                print("Note that using 'year' or 'year month' requires a column called 'date'.")
            print("\n")
        else:
            # Check to see if some categories are missing to warn the user
            group_by = set(['date' if cat in ['year', 'month'] else cat
                            for cat in self.discriminators])
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

            # subsample each groups, either by taking the sequences_per_group highest priority strains or
            # sampling at random from the sequences in the group
            seq_subsample = []
            for group, sequences_in_group in seq_names_by_group.items():
                if args.priority: #sort descending by priority
                    seq_subsample.extend(sorted(sequences_in_group, key=lambda x:priorities[x], reverse=True)[:spg])
                else:
                    seq_subsample.extend(sequences_in_group if len(sequences_in_group)<=self.sequences_per_group
                                         else random.sample(sequences_in_group, self.sequences_per_group))

            num_excluded_subsamp = len(self.sequence_names) - len(seq_subsample)
            return seq_subsample

    def divide_into_groups(self, sequences):
        groups = defaultdict(list)

        for sequence in sequences:
            group_name = self.group_name_for_sequence(sequence)
            groups[group_name].append(sequence)

        return groups

    def group_name_for_sequence(self, sequence):
        """
        Group names are roughly a tuple of the specified discriminators.
        "month" and "year" are special cases.
        """
        group_name = []
        for discriminator in self.discriminators:
            if discriminator in sequence.metadata:
                # if the discriminator is a metadata column, append it to the group name
                group_name.append(sequence.metadata[discriminator])
            elif discriminator in ["month", "year"] and "date" in sequence.metadata:
                # if the discriminator is "month" or "year", append the actual value to the group name
                try:
                    year = int(sequence.metadata["date"].split("-")[0])
                except ValueError:
                    # TODO does the sequence actually get skipped like the message says? Is that desired?
                    print(f"WARNING: Tried to --group-by year, but sample has no valid year. Skipping it. {sequence}")
                    continue

                if discriminator == "month":
                    try:
                        month = int(sequence.metadata["date"].split("-")[1])
                    except ValueError:
                        month = random.randint(1, 12)
                    group_name.append((year, month))
                else:
                    group_name.append(year)
            else:
                # otherwise, this discriminator causes an "unknown"
                group_name.append("unknown")

        return tuple(group_name)

    @property
    @functools.lru_cache()
    def priority_scores(self):
        # TODO raise exception if some samples do not have a priority? User could inadvertantly leave some samples unprioritized.
        priorities = defaultdict(float)

        try:
            with open(self.args.priority) as pfile:
                for name, priority in (
                    line.strip().split() for line in pfile.readlines()
                ):
                    priorities[name]: float(priority)
        except Exception as e:
            print(
                f"ERROR: missing or malformed priority scores file {self.args.priority}",
                file=sys.stderr,
            )
            raise e

        return priorities

    @classmethod
    def build(cls, context):
        return cls(
            priorities_fname=context.args.priority,
            group_by=context.args.group_by,
            sequences_per_group=context.args.sequences_per_group,
            seed=context.args.subsample_seed,
        )
