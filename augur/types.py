import enum


@enum.unique
class DataErrorMethod(enum.Enum):
    """
    Enum representation of string values that represent how a data error should
    be handled.
    """
    ERROR_FIRST     = 'error_first'
    ERROR_ALL       = 'error_all'
    WARN            = 'warn'
    SILENT          = 'silent'


@enum.unique
class ValidationMode(enum.Enum):
    """
    Enum representation of string values that represent how validation should
    be handled.
    """
    ERROR = 'error'
    WARN  = 'warn'
    SKIP  = 'skip'

    def __str__(self) -> str:
        """
        Stringify to the enum member's :py:attr:`.value` instead of the default.

        This let us use the enum's constructor and members with argparse's
        ``type`` and ``choices`` parameters, respectively, without exposing the
        enum class name to users.
        """
        return self.value
