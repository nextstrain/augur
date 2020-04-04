import factory
import random

import augur

import pytest


class SequenceFactory(factory.Factory):
    class Meta:
        model = augur.sequence.Sequence

    name = factory.Sequence(lambda n: f"seq{n}")
    sequence = str.join("", random.choices("ACGT", k=50))
    metadata = {}  # TODO
    bio_seq = factory.Sequence(lambda n: f"bio_seq{n}")


class SamplefileFactory(factory.Factory):
    class Meta:
        model = augur.sequence_file.SequenceFile
        exclude = ("sample_names",)

    filename = factory.Sequence(lambda n: f"seq{n}")
    metadata_filename = factory.Sequence(lambda n: f"seq{n}")


class VCFSelectorFactory(factory.Factory):
    class Meta:
        model = augur.filtering.vcf_selector.VCFSelector

    selected_samples = set(SequenceFactory.build_batch(3))
    samplefile = SamplefileFactory.create(filename="input.vcf")
    output_fname = "output.vcf"
