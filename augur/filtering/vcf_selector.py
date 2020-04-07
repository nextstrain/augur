import functools
import os

from augur.utils import run_shell_command, shquote


class VCFSelector:
    """
    Use vcftools to select the specified samples from `samplefile` and write an output file.
    """

    def __init__(self, *, selected_samples, samplefile, output_fname):
        self.selected_samples = selected_samples
        self.samplefile = samplefile
        self.output_fname = output_fname

    def write_output_vcf(self):
        print("Filtering samples using VCFTools with the call:")
        print(" ".join(self.vcftools_command))
        run_shell_command(" ".join(self.vcftools_command), raise_errors=True)
        # remove vcftools log file
        try:
            os.remove("out.log")
        except OSError:
            pass

    @property
    @functools.lru_cache()
    def vcftools_command(self):
        return (
            ["vcftools"]
            + self.vcftools_args_select_samples
            + self.vcftools_args_input_file
            + ["--recode", "--stdout"]
            + self.vcftools_args_output_redirect
        )

    @property
    def vcftools_args_select_samples(self):
        return [
            f(x)
            for x in self.selected_samples
            for f in (lambda _: "--indv", lambda sample: sample.name)
        ]

    @property
    def vcftools_args_input_file(self):
        fname = shquote(self.samplefile.filename)

        if self.samplefile.is_compressed_vcf:
            return ["--gzvcf", fname]
        else:
            return ["--vcf", fname]

    @property
    def vcftools_args_output_redirect(self):
        fname = shquote(self.output_fname)

        if fname.lower().endswith(".gz"):
            gzip_pipe = "| gzip -c"
        else:
            gzip_pipe = ""

        return [gzip_pipe, ">", fname]
