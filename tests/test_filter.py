import augur.filter


class TestFilter:
    def test_read_vcf_compressed(self):
        seq_keep, all_seq = augur.filter.read_vcf(
            "tests/builds/tb/data/lee_2015.vcf.gz"
        )

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq

    def test_read_vcf_uncompressed(self):
        seq_keep, all_seq = augur.filter.read_vcf("tests/builds/tb/data/lee_2015.vcf")

        assert len(seq_keep) == 150
        assert seq_keep[149] == "G22733"
        assert seq_keep == all_seq
