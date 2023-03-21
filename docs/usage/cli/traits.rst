============
augur traits
============

.. contents:: Table of Contents
   :local:

----

.. argparse::
    :module: augur
    :func: make_parser
    :prog: augur
    :path: traits
        


What about missing data?
========================

If you have strains with missing data and you want them to be reconstructed, then you must give them the value ``?``.
For example, if you are running a reconstruction of ``country`` and you don't know the country for a particular strain, you should set country to ``?`` in the metadata file for that strain.
Then, ``traits`` will estimate the most likely ``country`` value for any strains where you have provided ``?``. 

If you do not want these traits to be reconstructed (you would like it to remain clear that the ``country`` is unknown for this sample), then simply leave this field blank in the metadata file.

Note that each value -- empty strings, ``NA``, ``unknown``-- will be interpretted as a valid value!
So, it's best to be consistant with whatever you use for missing values, or those with ``NA`` will be shown as different from those with ``unknown``!
