=====================
Environment variables
=====================

Augur's behaviour can be globally modified by the values of some specific environment variables.
These can be especially useful in the context of an entire pipeline or workflow which uses Augur, as the environment variables can be set once for all Augur commands at the start of the pipeline.

.. envvar:: AUGUR_DEBUG

    Boolean.
    If set to a non-empty value, more detailed debugging information is shown by Augur during execution and handling of errors.

    For example, when this is enabled, stack traces and parent exceptions in an exception chain are no longer omitted from handled (i.e. expected) errors.
    Some commands will also emit more verbose operation logging during their execution.

.. envvar:: AUGUR_MINIFY_JSON

    Boolean.
    If set to a non-empty value, all JSON output produced by Augur will be minified by omitting indentation and newlines.

    Minifying the JSON will substantially reduce file sizes, which is helpful for large, deeply nested trees.

.. envvar:: AUGUR_RECURSION_LIMIT

    Integer.
    If set to a non-empty value, the Python recursion limit will be set to the given value early in Augur's execution by calling :func:`sys.setrecursionlimit`.

    Generally there is no need to set this environment variable.
    You may need to if you find yourself encountering :class:`RecursionError` while processing a very unbalanced tree.
