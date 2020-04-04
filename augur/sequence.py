import functools


# TODO move this
NUCLEOTIDE_SYMBOLS = {
    "A",
    "C",
    "G",
    "T",
    "-",
    "N",
    "R",
    "Y",
    "S",
    "W",
    "K",
    "M",
    "D",
    "H",
    "B",
    "V",
    "?",
}


# TODO rename this class to Sample
class Sequence:
    """
    A sequence sample.

    name : str
    sequence : str
    metadata : dict
        extra data not included in a FASTA or VCF file
    bio_seq : Bio.Seq.Seq
        optional original record from Bio, to be written out after filtering
    """

    def __init__(self, *, name, sequence, metadata, bio_seq=None):
        self.name = name
        self.sequence = sequence
        self.metadata = metadata
        self.bio_seq = bio_seq

    @property
    def date(self):
        return self.metadata.get("date")

    @property
    @functools.lru_cache()
    def length(self):
        return sum(
            1 for s in self.sequence if s in {"a", "t", "g", "c", "A", "T", "G", "C"}
        )

    @property
    @functools.lru_cache()
    def has_only_nucleotide_symbols(self):
        return any(site not in NUCLEOTIDE_SYMBOLS for site in self.sequence.upper())
