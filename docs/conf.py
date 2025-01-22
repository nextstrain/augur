# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('..'))


# -- Project information -----------------------------------------------------

from augur.__version__ import __version__ as augur_version
from datetime import date
import subprocess

def git_authors():
    result = subprocess.run(
        ["git", "shortlog", "--summary", "HEAD"],
        stdout = subprocess.PIPE,
        check  = True)

    names = [
        line.strip().split("\t")[1]
            for line in result.stdout.decode("utf-8").splitlines()
    ]

    return names

def prose_list(items):
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    elif len(items) == 2:
        return " and ".join(items)
    else:
        return ", ".join([*items[0:-1], "and " + items[-1]])

project = 'Augur'
version = augur_version
release = version
copyright = '2014–%d Trevor Bedford and Richard Neher' % (date.today().year)
author = prose_list(git_authors())


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'recommonmark',
    'sphinx.ext.autodoc',
    'sphinxarg.ext',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints', # must come after napoleon https://github.com/tox-dev/sphinx-autodoc-typehints/blob/1.21.4/README.md#compatibility-with-sphinxextnapoleon
    'sphinx_markdown_tables',
    'sphinx.ext.intersphinx',
    'sphinx_tabs.tabs',
    'nextstrain.sphinx.theme',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store',
    'contribute/DEV_DOCS.md',
    'faq/colors.md',
    'faq/fasta_input.md',
    'faq/import-beast.md',
    'faq/lat_longs.md',
    'faq/seq_traits.md',
    'faq/translate_ref.md',
    'faq/vcf_input.md',
    'usage/augur_snakemake.md',
]

# A string of reStructuredText that will be included at the end of every source
# file that is read. This is a possible place to add substitutions that should
# be available in every file.
rst_epilog = f"""
.. |authors| replace:: {author}
"""


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'nextstrain-sphinx-theme'

html_theme_options = {
    'logo_only': False, # if True, don't display project name at top of the sidebar
    'collapse_navigation': False, # if True, no [+] icons in sidebar
    'titles_only': True, # if True, page subheadings not included in nav
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# These paths are either relative to html_static_path
# or fully qualified paths (eg. https://...)
html_css_files = [
    'css/custom.css',
]

# -- Resolve build warnings --------------------------------------------------

nitpick_ignore = [
    # These are valid numpydoc keywords¹, but somehow they are not recognized by
    # napoleon.
    # ¹ https://numpydoc.readthedocs.io/en/v1.5.0/format.html#parameters
    ('py:class', 'optional'),
    ('py:class', 'iterable'),

     # Some references get translated to these, but somehow they can't get
     # resolved by intersphinx for a proper link.
     ("py:class", "json.decoder.JSONDecodeError"),
     ("py:class", "json.encoder.JSONEncoder"),

     # This class can't be referenced.
     # <https://github.com/python/cpython/issues/101503>
     ("py:class", "argparse._SubParsersAction"),
]

# -- Cross-project references ------------------------------------------------

intersphinx_mapping = {
    'Bio': ('https://biopython.org/docs/latest/', None),
    'docs.nextstrain.org': ('https://docs.nextstrain.org/en/latest/', None),
    'cli': ('https://docs.nextstrain.org/projects/cli/en/stable', None),
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
    'treetime': ('https://treetime.readthedocs.io/en/stable/', None),
}

# -- Linkchecking ------------------------------------------------------------

## NOTE: for both sets of regular expressions that follow, the
## underlying linkchecker code uses `re.match()` to apply them to URLs
## — so there's already an implicit "only at the beginning of a
## string" matching happening, and something like a plain `r'google'`
## regular expression will _NOT_ match all google.com URLs.
linkcheck_ignore = [
     # These URLs will occasionally fail and return 403 (broken).
     # <https://github.com/nextstrain/.github/issues/106#issuecomment-2408239782>
     r'^http://www\.microbesonline\.org/fasttree/',
     r'^https://academic\.oup\.com/ve/article/4/1/vex042/4794731',
     r'https://www\.gnu\.org/software/bash/manual/bash\.html#ANSI_002dC-Quoting',
]
linkcheck_anchors_ignore_for_url = [
     # Github uses anchor-looking links for highlighting lines but
     # handles the actual resolution with Javascript, so skip anchor
     # checks for Github URLs:
     r'^https://github\.com',
]
