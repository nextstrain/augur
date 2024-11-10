=============================
Examples of Augur in the wild
=============================

We really enjoy hearing about people using Augur.  If you use Augur to do your
work, we'd love to add it to this page!  Let us know about it by `opening an
issue <https://github.com/nextstrain/augur/issues/new?title=Augur%20in%20the%20wild>`__
or `submitting a PR <https://github.com/nextstrain/augur/pulls>`__.


Nextstrain pathogens
====================

All pathogens on `Nextstrain <https://nextstrain.org>`__ currently use Augur to
go from raw data to `Auspice <https://github.com/nextstrain/auspice>`__-ready
results files by performing cleanup, subsampling, alignment, tree building, and
ancestral state reconstruction.

Each pathogen has its own build repository on `GitHub <https://github.com>`__
which uses `Snakemake <http://snakemake.readthedocs.io>`__ to define and run
the pipeline of Augur commands.

=============================================================   =======================================================================================
Pathogen  (link to Nextstrain)                                  GitHub repository
=============================================================   =======================================================================================
`Seasonal influenza <https://nextstrain.org/flu/seasonal>`__    `nextstrain/seasonal-flu <https://github.com/nextstrain/seasonal-flu>`__
`West Nile Virus <https://nextstrain.org/WNV>`__                `nextstrain/WNV <https://github.com/nextstrain/WNV>`__
`Lassa <https://nextstrain.org/lassa>`__                        `nextstrain/lassa <https://github.com/nextstrain/lassa>`__
`Mumps <https://nextstrain.org/mumps>`__                        `nextstrain/mumps <https://github.com/nextstrain/mumps>`__
`Zika <https://nextstrain.org/zika>`__                          `nextstrain/zika <https://github.com/nextstrain/zika>`__
`Ebola <https://nextstrain.org/ebola>`__                        `nextstrain/ebola <https://github.com/nextstrain/ebola>`__
`Dengue <https://nextstrain.org/dengue>`__                      `nextstrain/dengue <https://github.com/nextstrain/dengue>`__
`Avian influenza <https://nextstrain.org/flu/avian>`__          `nextstrain/avian-flu <https://github.com/nextstrain/avian-flu>`__
`Measles <https://nextstrain.org/measles>`__                    `nextstrain/measles <https://github.com/nextstrain/measles>`__
`Enterovirus D68 <https://nextstrain.org/enterovirus/d68>`__    `neherlab/enterovirus_nextstrain <https://github.com/neherlab/enterovirus_nextstrain>`__
`Tuberculosis <https://nextstrain.org/tb>`__                    `nextstrain/tb <https://github.com/nextstrain/tb>`__
=============================================================   =======================================================================================

You can download a copy of a build repository, inspect the Augur commands used,
and try running it locally yourself.  This can be useful to reproduce the
analysis or tweak it for your own needs and data.


Nextstrain community
====================

Nextstrain also provides access to analyses by independent groups `stored and accessed via public GitHub repos <https://nextstrain.org/docs/contributing/community-builds>`__, all of which currently use augur.
For the most up-to-date listing of these please see the `community section on the nextstrain site. <https://nextstrain.org/community>`__
