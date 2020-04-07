from augur.filtering.vcf_selector import VCFSelector

import pytest


class TestVCFSelector:
    def test_vcf_command(
        self, vcf_selector_factory, samplefile_factory, sequence_factory
    ):
        vcf_selector = vcf_selector_factory.build(
            selected_samples=[sequence_factory.build(name=name) for name in ["a", "b"]],
            samplefile=samplefile_factory.build(filename="input.vcf.gz"),
            output_fname="output.vcf",
        )

        assert vcf_selector.vcftools_command == [
            "vcftools",
            "--indv",
            "a",
            "--indv",
            "b",
            "--gzvcf",
            "input.vcf.gz",
            "--recode",
            "--stdout",
            "",
            ">",
            "output.vcf",
        ]

    def test_vcf_args_select_samples(self, vcf_selector_factory, sequence_factory):
        vcf_selector = vcf_selector_factory.build(
            selected_samples=[sequence_factory.build(name=name) for name in ["a", "b"]]
        )

        assert vcf_selector.vcftools_args_select_samples == [
            "--indv",
            "a",
            "--indv",
            "b",
        ]

    @pytest.mark.parametrize(
        "test_fname, expected_option",
        [("input.vcf", "--vcf"), ("input.vcf.gz", "--gzvcf")],
    )
    def test_vcf_args_input_file(
        self, vcf_selector_factory, samplefile_factory, test_fname, expected_option
    ):
        vcf_selector = vcf_selector_factory.build(
            samplefile=samplefile_factory.build(filename=test_fname)
        )

        assert vcf_selector.vcftools_args_input_file == [expected_option, test_fname]

    @pytest.mark.parametrize(
        "test_fname, expected_pipe",
        [("output.vcf", ""), ("output.vcf.gz", "| gzip -c")],
    )
    def test_vcf_args_output_redirect(
        self, vcf_selector_factory, test_fname, expected_pipe
    ):
        vcf_selector = vcf_selector_factory.build(output_fname=test_fname)

        assert vcf_selector.vcftools_args_output_redirect == [
            expected_pipe,
            ">",
            test_fname,
        ]
