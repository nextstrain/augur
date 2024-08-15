"""
Merge two or more metadata tables into one.

Tables must be given unique names to identify them in the output and are
merged in the order given.

Rows are joined by id (e.g. "strain" or "name" or other
--metadata-id-columns), and ids must be unique within an input table (i.e.
tables cannot contain duplicate ids).  All rows are output, even if they
appear in only a single table (i.e. a full outer join in SQL terms).

Columns are combined by name, either extending the combined table with a new
column or overwriting values in an existing column.  For columns appearing in
more than one table, non-empty values on the right hand side overwrite values
on the left hand side.  The first table's id column name is used as the output
id column name.

One generated column per input table is appended to the end of the output
table to identify the source of each row's data.  Column names are generated
as "__source_metadata_{NAME}" where "{NAME}" is the table name given to
--metadata.  Values in each column are 1 or 0 for present or absent in that
input table.

Metadata tables of arbitrary size can be handled, limited only by available
disk space.  Tables are not required to be entirely loadable into memory.  The
transient disk space required is approximately the sum of the uncompressed size
of the inputs.

SQLite is used behind the scenes to implement the merge, but, at least for now,
this should be considered an implementation detail that may change in the
future.  The SQLite 3 CLI, sqlite3, must be available.  If it's not on PATH (or
you want to use a version different from what's on PATH), set the SQLITE3
environment variable to path of the desired sqlite3 executable.
"""
import gettext
import os
import re
import subprocess
import sys
from functools import reduce
from itertools import starmap
from shlex import quote as shquote
from shutil import which
from tempfile import mkstemp
from textwrap import dedent
from typing import Iterable, Tuple, TypeVar

from augur.argparse_ import ExtendOverwriteDefault, SKIP_AUTO_DEFAULT_IN_HELP
from augur.errors import AugurError
from augur.io.metadata import DEFAULT_DELIMITERS, DEFAULT_ID_COLUMNS, Metadata
from augur.io.print import print_err, print_debug
from augur.utils import first_line


T = TypeVar('T')


# Use ngettext() without a message catalog for its singular/plural handling so
# we can make proper error messages.  gettext() (no "n") is conventionally
# aliased as "_", so alias ngettext() as "_n".
_n = gettext.NullTranslations().ngettext


class NamedMetadata(Metadata):
    name: str
    """User-provided descriptive name for this metadata file."""

    table_name: str
    """Generated SQLite table name for this metadata file, based on *name*."""

    def __init__(self, name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name
        self.table_name = f"metadata_{self.name}"

    def __repr__(self):
        return f"<NamedMetadata {self.name}={self.path}>"


def register_parser(parent_subparsers):
    parser = parent_subparsers.add_parser("merge", help=first_line(__doc__))

    input_group = parser.add_argument_group("inputs", "options related to input")
    input_group.add_argument("--metadata", nargs="+", action="extend", required=True, metavar="NAME=FILE", help="Required. Metadata table names and file paths. Names are arbitrary monikers used solely for referring to the associated input file in other arguments and in output column names. Paths must be to seekable files, not unseekable streams. Compressed files are supported." + SKIP_AUTO_DEFAULT_IN_HELP)

    input_group.add_argument("--metadata-id-columns", default=DEFAULT_ID_COLUMNS, nargs="+", action=ExtendOverwriteDefault, metavar="COLUMN", help=f"Possible metadata column names containing identifiers, considered in the order given. Columns will be considered for all metadata tables. Only one ID column will be inferred for each table. (default: {' '.join(map(shquote_humanized, DEFAULT_ID_COLUMNS))})" + SKIP_AUTO_DEFAULT_IN_HELP)
    input_group.add_argument("--metadata-delimiters", default=DEFAULT_DELIMITERS, nargs="+", action=ExtendOverwriteDefault, metavar="CHARACTER", help=f"Possible field delimiters to use for reading metadata tables, considered in the order given. Delimiters will be considered for all metadata tables. Only one delimiter will be inferred for each table. (default: {' '.join(map(shquote_humanized, DEFAULT_DELIMITERS))})" + SKIP_AUTO_DEFAULT_IN_HELP)

    output_group = parser.add_argument_group("outputs", "options related to output")
    output_group.add_argument('--output-metadata', required=True, metavar="FILE", help="Required. Merged metadata as TSV. Compressed files are supported." + SKIP_AUTO_DEFAULT_IN_HELP)
    output_group.add_argument('--quiet', action="store_true", default=False, help="Suppress informational and warning messages normally written to stderr. (default: disabled)" + SKIP_AUTO_DEFAULT_IN_HELP)

    return parser


def run(args):
    print_info = print_err if not args.quiet else lambda *_: None

    # Parse --metadata arguments
    if not len(args.metadata) >= 2:
        raise AugurError(f"At least two metadata inputs are required for merging.")

    if unnamed := [repr(x) for x in args.metadata if "=" not in x or x.startswith("=")]:
        raise AugurError(dedent(f"""\
            All metadata inputs must be assigned a name, e.g. with NAME=FILE.

            The following {_n("input was", "inputs were", len(unnamed))} missing a name:

              {indented_list(unnamed, '            ' + '  ')}
            """))

    metadata = [name_path.split("=", 1) for name_path in args.metadata]

    if duplicate_names := [repr(name) for name, count
                                       in count_unique(name for name, _ in metadata)
                                       if count > 1]:
        raise AugurError(dedent(f"""\
            Metadata input names must be unique.

            The following {_n("name was", "names were", len(duplicate_names))} used more than once:

              {indented_list(duplicate_names, '            ' + '  ')}
            """))


    # Infer delimiters and id columns
    metadata = [
        NamedMetadata(name, path, args.metadata_delimiters, args.metadata_id_columns)
            for name, path in metadata]


    # Locate how to re-invoke ourselves (_this_ specific Augur).
    if sys.executable:
        augur = f"{shquote(sys.executable)} -m augur"
    else:
        # A bit unusual we don't know our own Python executable, but assume we
        # can access ourselves as the ``augur`` command.
        augur = f"augur"


    # Work with a temporary, on-disk SQLite database under a name we control so
    # we can access it from multiple (serial) processes.
    db_fd, db_path = mkstemp(prefix="augur-merge-", suffix=".sqlite")
    os.close(db_fd)

    # Clean up database file by default
    delete_db = True

    # Track columns as we see them, in order.  The first metadata's id column
    # is always the first output column of the merge, so insert it now.
    output_id_column = metadata[0].id_column
    output_columns = { output_id_column: [] }

    try:
        # Read all metadata files into a SQLite db
        for m in metadata:
            # All other metadata reading in Augur (i.e. via the csv module)
            # uses Python's "universal newlines"¹ definition and accepts \n,
            # \r\n, and \r as newlines interchangably (even mixed within the
            # same file!).  We accomplish the same behaviour here with SQLite's
            # less flexible newline handling by relying on the universal
            # newline translation of `augur read-file`.
            #   -trs, 24 July 2024
            #
            # ¹ <https://docs.python.org/3/glossary.html#term-universal-newlines>
            newline = os.linesep

            print_info(f"Reading {m.name!r} metadata from {m.path!r}…")
            sqlite3(db_path,
                f'.mode csv',
                f'.separator {sqlite_quote_dot(m.delimiter)} {sqlite_quote_dot(newline)}',
                f'.import {sqlite_quote_dot(f"|{augur} read-file {shquote(m.path)}")} {sqlite_quote_dot(m.table_name)}',

                f'create unique index {sqlite_quote_id(f"{m.table_name}_id")} on {sqlite_quote_id(m.table_name)}({sqlite_quote_id(m.id_column)});',

                # <https://sqlite.org/pragma.html#pragma_optimize>
                f'pragma optimize;')

            # We're going to use Metadata.columns to generate the select
            # statement, so ensure it matches what SQLite's .import created.
            assert m.columns == (table_columns := sqlite3_table_columns(db_path, m.table_name)), \
                f"{m.columns!r} == {table_columns!r}"

            # Track which columns appear in which metadata inputs, preserving
            # the order of both.
            for column in m.columns:
                # Match different id column names in different metadata files
                # since they're logically equivalent.
                output_column = output_id_column if column == m.id_column else column

                output_columns.setdefault(output_column, [])
                output_columns[output_column] += [(m.table_name, column)]


        # Construct query to produce merged metadata.
        select_list = [
            # Output metadata columns coalesced across input metadata columns
            *(f"""coalesce({', '.join(f"nullif({x}, '')" for x in starmap(sqlite_quote_id, reversed(input_columns)))}, null) as {sqlite_quote_id(output_column)}"""
                for output_column, input_columns in output_columns.items()),

            # Source columns
            *(f"""{sqlite_quote_id(m.table_name, m.id_column)} is not null as {sqlite_quote_id(f'__source_metadata_{m.name}')}"""
                for m in metadata)]

        from_list = [
            sqlite_quote_id(metadata[0].table_name),
            *(f"full outer join {sqlite_quote_id(m.table_name)} on {sqlite_quote_id(m.table_name, m.id_column)} in ({', '.join(sqlite_quote_id(m.table_name, m.id_column) for m in reversed(preceding))})"
                for m, preceding in [(m, metadata[:i]) for i, m in enumerate(metadata[1:], 1)])]

        # Take some small pains to make the query readable since it makes
        # debugging and development easier.  Note that backslashes aren't
        # allowed inside f-string expressions, hence the *newline* variable.
        newline = '\n'
        query = dedent(f"""\
            select
                {(',' + newline + '                ').join(select_list)}
            from
                {(newline + '                ').join(from_list)}
            ;
            """)


        # Write merged metadata as export from SQLite db.
        #
        # Assume TSV like nearly all other extant --output-metadata options.
        print_info(f"Merging metadata and writing to {args.output_metadata!r}…")
        print_debug(query)
        sqlite3(db_path,
            f'.mode csv',
            f'.separator "\\t" "\\n"',
            f'.headers on',
            f'.once {sqlite_quote_dot(f"|{augur} write-file {shquote(args.output_metadata)}")}',
            query)

    except SQLiteError as err:
        delete_db = False
        raise AugurError(str(err)) from err

    finally:
        if delete_db:
            os.unlink(db_path)
        else:
            print_info(f"WARNING: Skipped deletion of {db_path} due to error, but you may want to clean it up yourself (e.g. if it's large).")


def sqlite3(*args, **kwargs):
    """
    Internal helper for invoking ``sqlite3``, the SQLite CLI program.
    """
    sqlite3 = os.environ.get("SQLITE3", which("sqlite3"))

    if not sqlite3:
        raise AugurError(dedent(f"""\
            Unable to find the program `sqlite3`.  Is it installed?

            In order to use `augur merge`, the SQLite 3 CLI (version ≥3.39)
            must be installed separately.  It is typically provided by a
            Nextstrain runtime.
            """))

    argv = [sqlite3, "-batch", *args]

    print_debug(f"running {argv!r}")
    proc = subprocess.run(argv, encoding="utf-8", text=True, **kwargs)

    try:
        proc.check_returncode()
    except subprocess.CalledProcessError as err:
        raise SQLiteError(f"sqlite3 invocation failed") from err

    return proc


class SQLiteError(Exception):
    pass


def sqlite3_table_columns(db_path, table: str) -> Iterable[str]:
    return sqlite3(db_path, f"select name from pragma_table_info({sqlite_quote_string(table)})", capture_output=True).stdout.splitlines();


def sqlite_quote_id(*xs):
    """
    Quote a SQLite identifier.

    <https://sqlite.org/lang_keywords.html>

    >>> sqlite_quote_id('foo bar')
    '"foo bar"'
    >>> sqlite_quote_id('table name', 'column name')
    '"table name"."column name"'
    >>> sqlite_quote_id('weird"name')
    '"weird""name"'
    """
    return '.'.join('"' + x.replace('"', '""') + '"' for x in xs)


def sqlite_quote_dot(x):
    """
    Quote a SQLite CLI dot-command argument.

    <https://sqlite.org/cli.html#dot_command_arguments>
    """
    return '"' + x.replace('\\', '\\\\').replace('"', '\\"') + '"'


def sqlite_quote_string(x):
    """
    Quote a SQLite string (i.e. produce a string literal).

    <https://www.sqlite.org/lang_expr.html#literal_values_constants_>
    """
    return "'" + x.replace("'", "''") + "'"


def count_unique(xs: Iterable[T]) -> Iterable[Tuple[T, int]]:
    # Using reduce() with a dict because it preserves input order, unlike
    # itertools.groupby(), which requires a sort.  Preserving order is a nice
    # property for the user since we generate an error message with this.
    #   -trs, 24 July 2024
    yield from reduce(lambda counts, x: {**counts, x: counts.get(x, 0) + 1}, xs, counts := {}).items() # type: ignore


def indented_list(xs, prefix):
    return f"\n{prefix}".join(xs)


def shquote_humanized(x):
    r"""
    shquote for humans.

    Use C-style escapes supported by shells (specifically, Bash) for characters
    that humans would typically use C-style escapes for instead of quoted
    literals.

    <https://www.gnu.org/software/bash/manual/bash.html#ANSI_002dC-Quoting>

    >>> shquote_humanized("abc")
    'abc'

    >>> shquote_humanized("\t")
    "$'\\t'"

    >>> shquote_humanized("abc def")
    "'abc def'"

    >>> shquote_humanized("abc\tdef")
    "abc$'\\t'def"
    """
    escapes = {
        '\a': r'\a',
        '\b': r'\b',
        '\f': r'\f',
        '\n': r'\n',
        '\r': r'\r',
        '\t': r'\t',
        '\v': r'\v',
    }

    def quote(s):
        if s in escapes:
            return f"$'{escapes[s]}'"
        else:
            # split leaves leading and trailing empty strings when its input is
            # entirely (captured) separator.  Avoid quoting every empty string
            # *part* here…
            return shquote(s) if s else ''

    parts = re.split('([' + ''.join(escapes.values()) + '])', x)
    quoted = ''.join(map(quote, parts))

    # …and instead quote a final empty string down here if we're still empty
    # after joining all our parts together.
    return quoted if quoted else shquote('')
