============
augur traits
============

.. contents::

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: traits
        


What about missing data?
========================

If you have strains with missing data then you must give them the value ``?``.
For example, if you are running a reconstruction of ``country`` and you don't know the country for a particular strain, you must set country to ``?`` in the metadata file for that strain.
Note that anything else -- empty strings, ``NA``, ``unknown``-- will be interpretted as a valid value!

Currently there is no way to *not* infer the value of these missing data, but we are working on making this option available in the future.
