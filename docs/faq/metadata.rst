How do I prepare metadata?
==========================

Analyses are vastly more interesting if the sequences or samples
analyzed have rich 'meta data' wherever possible. This metadata could
typically include collection dates, geographic location, symptoms of
patients, host characteristics, etc.

To make the most of augur's features, we recommend including sampling
date and at least one type of geographic information if at all possible.
However, you can also include things like symptoms, host, clinical
outcome - and more!

For augur to be able to parse this data, it needs to be formated
consistently. Your data may have meta information coded into the
sequence name (see example :ref:`below<parsing-from-the-header>`). If
not, a very transparent way is to provide the meta data as a separate
table in a tab- or comma-separated file.

An example meta data file is shown here:

::

   strain      accession   date        region          host
   1_0087_PF   KX447509    2013-12-XX  Oceania         Human
   1_0181_PF   KX447512    2013-12-XX  Oceania         Bat
   1_0199_PF   KX447519    2013-11-XX  Oceania         Human
   BRA/2016    KY785433    2016-04-08  South America   Cow
   BRA/2015    KY558989    2015-02-23  South America   Bat

A note on Excel
~~~~~~~~~~~~~~~

Because Excel will automatically change the date formatting, we
recommend *not* opening or preparing your meta data file in Excel. If
the metadata is already in Excel, or you decide to prepare it in Excel,
we recommend using another program to correct the dates afterwards (and
don't open it in Excel again!).

Format
~~~~~~

**Strain names**

You must have one column named ``strain`` or ``name``. It contains your
sequence names, and needs to match the identifiers of your sequences (in
the Fasta or VCF file) *exactly* and must not contain characters such as
spaces, or ``()[]{}|#><``.

**Dates**

Dates should be formated according as ``YYYY-MM-DD``. You can specify
unknown dates or month by replacing the respected values by ``XX`` (ex:
``2013-01-XX`` or ``2011-XX-XX``) and completely unknown dates can be
shown with ``XXXX-XX-XX``.

**Geography**

Geographic locations can be broken down, for example, into ``region``,
``country``, ``division`` or ``city``. You can have as many levels of
geographic information as you wish. For ``region``, ``country``, and
some ``division``\ s augur already knows many lat-long coordinates (see
which ones it already knows by checking the list
`here <https://github.com/nextstrain/augur/blob/-/augur/data/lat_longs.tsv>`__).

It is important that these are spelled consistently.

If you want to include locations where augur doesn't know the lat-long
values, you can include them - see how :doc:`here <docs.nextstrain.org:guides/bioinformatics/lat_longs>`.

Consistancy and Style
~~~~~~~~~~~~~~~~~~~~~

Check that your metadata is free from spelling mistakes and that values
are consistant. Augur doesn't know that 'UK' and 'United Kingdom' or
'cat' and 'feline' are the same!

Previously, auspice 'prettified' traits by capitalizing them
automatically, and removing the underscores that separated two-word
locations ('new_zealand' became 'New Zealand').

Auspice will still do this if you are exporting 'V1' type JSON files
(from augur v5 or augur v6 using ``export v1``), but will not do this if
you are using ``export v2`` (:ref:`read more <prettifying-metadata-fields>`).
Instead, you should update your metadata files so that traits look the
same as you'd like them to display in Auspice (change 'new_zealand' to
'New Zealand' in your metadata, and in any additional latitude-longitude
or coloring files you use).

.. _parsing-from-the-header:

Parsing from the header
~~~~~~~~~~~~~~~~~~~~~~~

Sometimes, metadata can be coded into the Fasta header, like so:

::

   >1_0087_PF | KX447509 | 2013-12-XX | oceania
   ACTCGCTGCATCG...

Augur can parse meta data from Fasta headers using the ``parse``
function (see :doc:`here </usage/cli/parse>`), but you have to make sure
that every sequence has the exact same meta data fields (even if empty),
and that they are consistently delimited with ``|``. Furthermore, none
of the metadata fields can contain the character ``|``.
