from argparse import ArgumentTypeError
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

    @classmethod
    def argtype(cls, input_string):
        """
        Intended to be used as the argument type converter for argparse options
        that use the enum values as inputs.

        Raises a custom `argparse.ArgumentTypeError` so that the error
        message can include a helpful list of the all valid enum values.
        """
        try:
            return cls(input_string)
        except ValueError as error:
            choices = ", ".join(f"{str(x)!r}" for x in cls)
            raise ArgumentTypeError(
                f"invalid choice: {input_string!r} (choose from {choices})") from error


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
