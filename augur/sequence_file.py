from Bio import SeqIO
import functools

from augur.sequence import Sequence
from augur.utils import read_metadata


# TODO rename this to SampleFile
class SequenceFile:
    """
    Represents a file containing sequences (TODO or samples?) in FASTA, VCF, or gzipped VCF format.
    """
    def __init__(self, filename, metadata_filename):
        self.filename = filename
        self.metadata_filename = metadata_filename

    @property
    @functools.lru_cache()
    def sequences(self):
        """
        An array of Sequence objects representing the input sequences adorned with input metadata.
        Input sequences that lack metadata are not included.
        """
        if self.is_vcf:
            # Note that VCFs do not contain actual sequences, so we'll use None
            return [
                Sequence(name=name, sequence=None, metadata=self.metadata[name])
                for name in self.read_sequence_names_from_vcf()
                if name in self.metadata
            ]
        elif self.is_fasta:
            return [
                Sequence(
                    name=bio_seq.id,
                    sequence=bio_seq.seq,
                    metadata=self.metadata[bio_seq.id],
                    bio_seq=bio_seq,
                )
                for bio_seq in SeqIO.parse(self.filename, "fasta")
                if bio_seq.id in self.metadata
            ]
        else:
            raise  # TODO specific exception class for unsupported format

    def read_sequence_names_from_vcf(self):
        with self.open_vcf() as file:
            chrom_line = next(line for line in file if line.startswith("#C"))
            headers = chrom_line.strip().split("\t")

            return headers[headers.index("FORMAT") + 1 :]

    def open_vcf(self):
        if self.is_compressed_vcf:
            import gzip

            return gzip.open(self.filename, mode="rt")

        return open(self.filename)

    @property
    @functools.lru_cache()
    def metadata(self):
        meta_dict, _meta_columns = read_metadata(self.metadata_filename)
        return meta_dict

    @property
    def is_vcf(self):
        return self.is_compressed_vcf or self.is_uncompressed_vcf

    @property
    def is_compressed_vcf(self):
        return self.lower_filename.endswith(".vcf.gz")

    @property
    def is_uncompressed_vcf(self):
        return self.lower_filename.endswith(".vcf")

    @property
    def is_fasta(self):
        # TODO is this reliable?
        return self.lower_filename.endswith(".fasta")

    @property
    def lower_filename(self):
        return self.filename.lower()
