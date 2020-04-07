class AugurException(Exception):
    pass


class MissingDependencyException(AugurException):
    """Raised when a required external dependency is not present."""
    #TODO raise "ERROR: `vcftools` is not installed! This is required for VCF data. Please see the augur install instructions to install it."
    pass
