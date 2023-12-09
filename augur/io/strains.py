from augur.utils import read_entries


def read_strains(*files, comment_char="#"):
    """Reads strain names from one or more plain text files and returns the
    set of distinct strains.

    Strain names can be commented with full-line or inline comments. For
    example, the following is a valid strain names file::

        # this is a comment at the top of the file
        strain1  # exclude strain1 because it isn't sequenced properly
        strain2
          # this is an empty line that will be ignored.

    Parameters
    ----------
    files : iterable of str
        one or more names of text files with one strain name per line

    Returns
    -------
    set :
        strain names from the given input files

    """
    return set(read_entries(*files, comment_char=comment_char))
