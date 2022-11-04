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
