import enum


class ArgparseEnum(enum.Enum):
    """
    Intended to be used as a parent class for any enum representation of
    string values to be used with argparse options.

    Can be replaced by :py:class:`enum.StrEnum` once Augur's minimum supported
    Python version is 3.11.
    """
    def __str__(self) -> str:
        """
        Stringify to the enum member's :py:attr:`.value` instead of the default.

        This let us use the enum's constructor and members with argparse's
        ``type`` and ``choices`` parameters, respectively, without exposing the
        enum class name to users.
        """
        return self.value


@enum.unique
class DataErrorMethod(ArgparseEnum):
    """
    Enum representation of string values that represent how a data error should
    be handled.
    """
    ERROR_FIRST     = 'error_first'
    ERROR_ALL       = 'error_all'
    WARN            = 'warn'
    SILENT          = 'silent'


@enum.unique
class EmptyOutputReportingMethod(ArgparseEnum):
    """
    Enum representation of string values that represent how empty outputs should
    be reported.
    """
    ERROR = 'error'
    WARN  = 'warn'
    SILENT  = 'silent'


@enum.unique
class ValidationMode(ArgparseEnum):
    """
    Enum representation of string values that represent how validation should
    be handled.
    """
    ERROR = 'error'
    WARN  = 'warn'
    SKIP  = 'skip'
