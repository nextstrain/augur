import argparse
import json

import pytest

from augur.errors import AugurError
from augur.io.nextclade_dataset import read_nextclade_dataset, apply_nextclade_dataset


@pytest.fixture
def dataset(tmp_path):
    """A Nextclade dataset with a reference and a genome annotation in a subdirectory."""
    (tmp_path / "reference.fasta").write_text(">ref\nACGT\n")
    (tmp_path / "annotation").mkdir()
    (tmp_path / "annotation" / "genome_annotation.gff3").write_text("##gff-version 3\n")
    write_pathogen_json(tmp_path, {"files": {"reference": "reference.fasta", "genomeAnnotation": "annotation/genome_annotation.gff3"}})
    return tmp_path


def write_pathogen_json(path, content):
    (path / "pathogen.json").write_text(json.dumps(content))


class TestReadNextcladeDataset:
    def test_directory(self, dataset):
        result = read_nextclade_dataset(str(dataset))
        assert result.reference == str(dataset / "reference.fasta")
        assert result.annotation == str(dataset / "annotation" / "genome_annotation.gff3")

    def test_pathogen_json(self, dataset):
        result = read_nextclade_dataset(str(dataset / "pathogen.json"))
        assert result.reference == str(dataset / "reference.fasta")

    def test_without_annotation(self, dataset):
        write_pathogen_json(dataset, {"files": {"reference": "reference.fasta"}})
        assert read_nextclade_dataset(str(dataset)).annotation is None

    def test_missing_pathogen_json(self, tmp_path):
        with pytest.raises(AugurError, match="does not contain a pathogen.json"):
            read_nextclade_dataset(str(tmp_path))

    def test_invalid_pathogen_json(self, tmp_path):
        (tmp_path / "pathogen.json").write_text("{")
        with pytest.raises(AugurError, match="Could not parse"):
            read_nextclade_dataset(str(tmp_path))

    def test_missing_files_section(self, dataset):
        write_pathogen_json(dataset, {"schemaVersion": "3.0.0"})
        with pytest.raises(AugurError, match="does not contain a 'files' section"):
            read_nextclade_dataset(str(dataset))

    def test_missing_reference(self, dataset):
        write_pathogen_json(dataset, {"files": {"genomeAnnotation": "annotation/genome_annotation.gff3"}})
        with pytest.raises(AugurError, match=r"does not declare a reference sequence \(files.reference\)"):
            read_nextclade_dataset(str(dataset))

    def test_declared_file_missing(self, dataset):
        write_pathogen_json(dataset, {"files": {"reference": "reference.fasta", "genomeAnnotation": "missing.gff3"}})
        with pytest.raises(AugurError, match="'missing.gff3' declared as files.genomeAnnotation .* does not exist"):
            read_nextclade_dataset(str(dataset))


class TestApplyNextcladeDataset:
    def test_no_dataset(self):
        args = argparse.Namespace(nextclade_dataset=None, annotation="explicit.gff")
        apply_nextclade_dataset(args, annotation="annotation")
        assert args.annotation == "explicit.gff"

    def test_fills_arguments(self, dataset):
        args = argparse.Namespace(nextclade_dataset=str(dataset), vcf_reference=None, annotation=None, use_nextclade_gff_style=False)
        apply_nextclade_dataset(args, reference="vcf_reference", annotation="annotation")
        assert args.vcf_reference == str(dataset / "reference.fasta")
        assert args.annotation == str(dataset / "annotation" / "genome_annotation.gff3")
        assert args.use_nextclade_gff_style is True

    def test_reference_only(self, dataset):
        args = argparse.Namespace(nextclade_dataset=str(dataset), reference_sequence=None)
        apply_nextclade_dataset(args, reference="reference_sequence")
        assert args.reference_sequence == str(dataset / "reference.fasta")
        assert not hasattr(args, "use_nextclade_gff_style")

    def test_conflicting_argument(self, dataset):
        args = argparse.Namespace(nextclade_dataset=str(dataset), vcf_reference="explicit.fasta")
        with pytest.raises(AugurError, match="--nextclade-dataset and --vcf-reference can not be used together"):
            apply_nextclade_dataset(args, reference="vcf_reference")

    def test_annotation_required(self, dataset):
        write_pathogen_json(dataset, {"files": {"reference": "reference.fasta"}})
        args = argparse.Namespace(nextclade_dataset=str(dataset), annotation=None)
        with pytest.raises(AugurError, match=r"does not declare a genome annotation \(files.genomeAnnotation\)"):
            apply_nextclade_dataset(args, annotation="annotation")
