from textwrap import dedent

import pytest
from Bio.SeqFeature import CompoundLocation

from augur.errors import AugurError
from augur.io.sequences import load_features


HEADER = "##gff-version 3\n##sequence-region seq 1 1000\n"


def write_gff(tmp_path, rows):
    """Write a GFF with the given rows (whitespace-separated columns) and return its path."""
    path = tmp_path / "annotation.gff"
    lines = ["\t".join(row.split()) for row in dedent(rows).strip().splitlines()]
    path.write_text(HEADER + "\n".join(lines) + "\n")
    return str(path)


def segments(feature):
    return [(int(part.start), int(part.end), part.strand) for part in feature.location.parts]


class TestLoadFeaturesNextclade:
    def test_nuc(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . gene 1 9 . + . gene=A
        """)
        features = load_features(gff, nextclade_gff=True)
        assert (int(features['nuc'].location.start), int(features['nuc'].location.end)) == (0, 1000)

    def test_multi_row_cds_joined_in_file_order(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . gene 1 100 . + . ID=gene-A;gene=A
            seq . mRNA 1 100 . + . ID=rna-A;Parent=gene-A
            seq . CDS  31 39 . + 0 ID=cds-X;Parent=rna-A;Name=X
            seq . CDS  1  9  . + 0 ID=cds-X;Parent=rna-A;Name=X
        """)
        features = load_features(gff, nextclade_gff=True)
        assert set(features) == {'nuc', 'X'}
        assert isinstance(features['X'].location, CompoundLocation)
        assert segments(features['X']) == [(30, 39, 1), (0, 9, 1)]
        assert features['X'].type == 'CDS'
        assert features['X'].id == 'cds-X'

    def test_multi_row_cds_minus_strand(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 31 39 . - 0 ID=cds-X;Name=X
            seq . CDS 1  9  . - 0 ID=cds-X;Name=X
        """)
        features = load_features(gff, nextclade_gff=True)
        assert segments(features['X']) == [(30, 39, -1), (0, 9, -1)]

    def test_multi_row_cds_mixed_strands(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 1  9  . + 0 ID=cds-X;Name=X
            seq . CDS 31 39 . - 0 ID=cds-X;Name=X
        """)
        with pytest.raises(AugurError, match="rows on different strands"):
            load_features(gff, nextclade_gff=True)

    def test_cds_rows_without_id_are_separate(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 1  9  . + 0 Name=X
            seq . CDS 31 39 . + 0 Name=Y
        """)
        features = load_features(gff, nextclade_gff=True)
        assert segments(features['X']) == [(0, 9, 1)]
        assert segments(features['Y']) == [(30, 39, 1)]

    def test_gene_fallback_is_per_gene(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . gene 1  30 . + . ID=gene-A;gene=A
            seq . CDS  1  30 . + 0 ID=cds-A;Parent=gene-A;Name=A-cds
            seq . gene 31 60 . + . ID=gene-B;gene=B
        """)
        features = load_features(gff, nextclade_gff=True)
        assert set(features) == {'nuc', 'A-cds', 'B'}
        assert features['B'].type == 'gene'

    def test_name_priority(self, tmp_path):
        # CDSs prefer 'Name' over 'gene', genes prefer 'gene' over 'Name'
        gff = write_gff(tmp_path, """
            seq . CDS  1  30 . + 0 ID=cds-1;Name=cds-name;gene=cds-gene
            seq . gene 31 60 . + . ID=gene-1;Name=gene-name;gene=gene-gene
            seq . CDS  61 90 . + 0 ID=cds-2;product=cds-product;protein_id=P1
            seq . CDS  91 99 . + 0 ID=cds-3
        """)
        features = load_features(gff, nextclade_gff=True)
        assert set(features) == {'nuc', 'cds-name', 'gene-gene', 'cds-product', 'cds-3'}

    def test_unnamed_feature(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 1 9 . + 0 Note=nameless
        """)
        with pytest.raises(AugurError, match="without any of the attributes used for naming"):
            load_features(gff, nextclade_gff=True)

    def test_duplicate_names(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS  1  9  . + 0 ID=cds-1;Name=X
            seq . gene 31 39 . + . ID=gene-1;gene=X
        """)
        with pytest.raises(AugurError, match="multiple genes/CDSs with the name 'X'"):
            load_features(gff, nextclade_gff=True)

    def test_nuc_name(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 1 9 . + 0 Name=nuc
        """)
        with pytest.raises(AugurError, match="with the name 'nuc'"):
            load_features(gff, nextclade_gff=True)

    def test_no_cds_or_genes(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . region 1 1000 . + . ID=seq
        """)
        with pytest.raises(AugurError, match="contains no CDS or gene features"):
            load_features(gff, nextclade_gff=True)

    def test_feature_names(self, tmp_path, capsys):
        gff = write_gff(tmp_path, """
            seq . CDS 1  9  . + 0 Name=X
            seq . CDS 31 39 . + 0 Name=Y
        """)
        features = load_features(gff, feature_names=['Y', 'Z'], nextclade_gff=True)
        assert set(features) == {'nuc', 'Y'}
        assert "Couldn't find gene/CDS Z" in capsys.readouterr().out

    def test_length_validated_across_segments(self, tmp_path):
        gff = write_gff(tmp_path, """
            seq . CDS 1  10 . + 0 ID=cds-X;Name=X
            seq . CDS 31 34 . + 0 ID=cds-X;Name=X
        """)
        with pytest.raises(AugurError, match="'X' has length 14 which is not a multiple of 3"):
            load_features(gff, nextclade_gff=True)
