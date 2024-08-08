=============================
Migrating from augur v5 to v6
=============================

Augur is a versatile bioinformatics toolkit in its own right, but it
is developed in conjunction with `Auspice
<https://docs.nextstrain.org/projects/auspice/en/latest/>`__, the
interactive visualisation tool behind what you see on `nextstrain.org
<https://nextstrain.org>`__. Because of this, build pipelines often
finish with ``augur export`` which produces JSONs that Auspice can
visualise.

As Auspice grows and improves, the format of the JSON files is
changing too, to allow more functionality and flexibility. This means
Augur's export function needs to change as well. This page details how
the new ``augur export v2`` works and what you need to do to start
using it.

Some important points:

* 'Old' JSON files (v1, made with Augur v5) still work with the new
  (v2) Auspice for the time being - but may not be supported in future
  versions.
* When you upgrade to Augur v6, ``augur export`` will no longer work.
  You'll need to specify which version of Auspice (v1 or v2) your JSON
  files should be.

----

.. contents::
   :local:
   :depth: 1

----

Compatibility between Augur & Auspice versions
----------------------------------------------

Augur v5 simply used ``augur export`` to produce two JSONs (tree +
meta) for Auspice v1 - so we'll call these JSONs "v1" JSONs.

The new Augur (v6) can still create "v1" JSONS, but can also create
JSONs that work with the latest `Auspice release
<https://docs.nextstrain.org/projects/auspice/en/latest/releases/v2.html>`__ -
Auspice v2. The new format combines the tree and metadata into one
JSON, which we'll call a "v2" JSON.

:doc:`This page <./auspice-compatibility>` has the most
up-to-date compatibility information between different Augur and
Auspice versions.

We understand how important backwards compatibility is - so for the
time being "v1" JSONs will continue to work with Auspice v2. However,
we recommend switching to v2 JSONs - they have more features, are
easier to work with, and future versions of Auspice may not support v1
JSONs!


Motivation behind changing JSON formats
---------------------------------------

With the release of `Auspice v2
<https://docs.nextstrain.org/projects/auspice/en/latest/releases/v2.html>`__,
we made some breaking changes to the JSON schema. Why change formats?
We were motivated by:

* **Compactness**: Tree and Meta JSON files are now combined, so you
  only have to worry about one output file
* **Flexibility**: The new v2 JSONs allow us flexibility to include
  more features and data, and will let us move towards getting in line
  with existing conventions like GFF and BibTex
* **Ease of use**: Users commonly got confused by the 'config' file.
  For basic runs you can now specify everything you need to see your
  data right in the command-line - no 'config' file needed! For more
  advanced exports, you can still specify a config file with more
  detail. (See ":ref:`Using a Config File <using-a-config-file-to-customise-the-visualisation>`")

I just need my old run to work *right now*!
-------------------------------------------

We completely understand you may not have time to make the change
right this second, and that there's nothing more frustrating than
having a run break right before a presentation or deadline!

If you want to keep using the old version *for now*, replace ``augur
export`` with ``augur export v1`` - everything else remains the same.
You can use Auspice v1 or Auspice v2 (see :doc:`compatibility <./auspice-compatibility>`).

To use the new version, use ``augur export v2``. You'll need to make a
few small changes, but you'll be future-proofing your runs. *(Future
you will thank past you!)*

What needs changing to use ``augur export v2``?
-----------------------------------------------

.. note:: Helpful hint: You can always get a full overview of the
          arguments for export v2 with ``augur export v2 --help``

What's the same
***************

You still pass in your tree, metadata, and node-data files with
``--tree``, ``--metadata``, and ``--node-data`` - just like in ``export v1``.
Similarly, you can pass in files containing colors, latitude, and
longitude data using ``--colors`` and ``--lat-longs``, respectively.

.. _different-outputs:

Different outputs
*****************

Instead of specifying two output files (``--output-tree`` and
``--output-meta``) you now only need to specify one with ``--output``.
For example, if your old files were ``auspice/virus_AB_tree.json`` and
``auspice/virus_AB_meta.json``, you might want to call the single
output ``auspice/virus_AB.json`` - or if you want to tell it apart
from your v1 export, you might call it ``auspice/virus_ABv2.json``.

To export the reference sequence relative to which mutations have been
identified, specify the ``--include-root-sequence`` flag. This flag
writes a JSON whose name is relative to the stem of the main output
JSON. For VCF input, this file will contain the reference sequence to
which the VCF is mapped. For example, if the main output is called
``auspice/virus_AB.json``, the root sequence will be saved to
``auspice/virus_AB_root-sequence.json``.

Other changed arguments
***********************

The ``--tree-name`` argument has been removed, as auspice v2 no longer
uses this. See `the auspice docs
<https://docs.nextstrain.org/projects/auspice/en/latest/advanced-functionality/second-trees.html>`__
for more information about how second trees are specified and
displayed.

Command Line Options instead of (or in addition to) a config file
*****************************************************************

One of the biggest changes in ``augur export v2`` is that you can pass
much more in using the command-line, meaning 'config' files are no
longer required. The 'config' or 'Auspice config' file defines a
number of visualisation settings such as title, default displays, and
which colorings to use. However, it's been a source of pain for many
users!

Many of these things can now be passed in by the command-line, but
some options are only possible using the config file. You can always
continue to put most things in the config file if you prefer. If you
want to use a :ref:`config file
<using-a-config-file-to-customise-the-visualisation>`, you can pass
this in with ``--auspice-config``, but the format of this has changed
(see the link).

It's important to note that generally *any command line options you
use will override the same option in your config file*.

Coloring traits is smarter
**************************

Previously, anything you wanted to color by had to be in the config
file. You always had to include a ``gt`` and ``num_date`` entry, and
remember to add anything new to the file.

We've made this smarter - ``augur export v2`` now automatically detects
some traits and you can specify others on the command line. You can
also control the color options in more detail using a config file.

We'll cover how coloring works on the :ref:`command line
<command-line>` and how it works in :ref:`config files <colorings>` in more detail
below.

Traits display exactly how you want
***********************************

Previously, auspice tried to make traits and locations look 'pretty'
by auto-capitalizing them and removing underscores (which were
required in multi-word traits). Auspice no longer does this for v2
JSONS, so you'll need to ensure your traits look exactly how you want
them to display in auspice. You can read more about that
:ref:`here <prettifying-metadata-fields>`.

-----------------------------------------------------------------

Terminology
-----------

.. _terminology-traits:

Traits
******

Traits is the general term for certain data associated with nodes in
the tree, for example "country", "serotype", or "age". These may have
been inferred for internal nodes by Augur functions like ``augur
traits`` (confusingly named!) and ``augur clades``, or they may only be
available for tips and provided by the metadata TSV file.

Geographic Traits
*****************

Certain traits have a geographic interpretation, e.g. "country".
Auspice will attempt to display these traits on a map (and provide a
drop-down to switch between them if there are more than one).

.. note:: Make sure that these have corresponding entry in the
          lat-longs TSV file supplied to ``export``. See how to do this
          `here <https://docs.nextstrain.org/en/latest/guides/bioinformatics/lat_longs.html>`_.

.. _command-line:

Using command-line options to customise the visualisation
---------------------------------------------------------

As mentioned above, you can now replace most of the functionality of
the "Auspice config" file with command line options. We hope that for
most users this means the config file isn't necessary (but it's always
there is you need its advanced functionality).

.. note:: Remember that generally *any command line options you use
          will override the same option in your config file*.

Title
*****

Set the title displayed by Auspice via ``--title`` (previously this
was the "title" field in your v1 config file). If running directly
from the command line, put your title in quotes (ex: ``--title
"Phylodynamics of my Pathogen"``). If you are using Snakemake and
passing the value using ``params``, you'll need to double-quote the
title using single and double quotes. For example:

.. code-block:: python

   params:
     title = "'Phylodynamics of my Pathogen'"
   shell:
     "augur export v2 --title {params.title} ..."

Maintainers
***********

The maintainer(s) are displayed in the footer of Auspice and may have
associated links. These can be specified with ``--maintainers`` and
you can now have more than one maintainer associated with your run.
Previously this was set by the "maintainer" field in your v1 config
file and was limited to a single entry.

If running directly from the command line, put each maintainer in
quotes (ex: ``--maintainers "Jane Doe" "Ravi Kupra"``). If you have a
URL associated with a maintainer (completely optional), then you can
add them like so:

.. code-block:: shell

   --maintainers "Jane Doe <mailto:jane.doe@example.com>" "Ravi Kupra <https://github.com/ravikupra"

If you are using Snakemake and passing the value using ``params``,
you'll need to put the whole list in double quotes, and each person in
single quotes. For example:

.. code-block:: python

   params:
     maints = "'Jane Doe' 'Ravi Kupra <github.com/ravikupra>'"
   shell:
     "augur export v2 --maintainers {params.maints} ..."

You will need to use quotes in the same way even if you only have one
maintainer!

Build URL
*********

Set the build URL displayed by Auspice via ``--build-url``. If running
directly from the command line, input your build URL directly (ex:
``--build-url https://github.com/nextstrain/zika``).

Description
***********

Set the description and/or acknowledgments in the footer of Auspice
via ``--description``. This option expects a Markdown file containing
description to be displayed.

Panels
******

Auspice will, by default, try to show the tree, map, and entropy
panels. You can customise this with the ``--panels`` option, which was
previously the "panels" field in the your v1 config file. Options are
"tree", "map", "entropy", and "frequencies" (e.g: ``--panels tree map
entropy``).

.. note:: If you want to display the frequencies panel, you must
          specify "frequencies" *and* ensure a tip frequency file is
          available for ``auspice`` to access.

.. _cli-traits:

Traits
******

(:ref:`What's a trait? <terminology-traits>`) Traits will become coloring options
in Auspice. Some are automatically included, and some can be defined
on the command line. The following rules are followed for which traits
will be exported:

1. Genotype and date (if present) are always automatically included as
   coloring options - you don't need to include them. *(Previously
   these were "gt" and "numdate" in the "color_options" section of
   your v1 config file.)*
2. Traits contained in the node-data JSONs handed to ``augur export``
   (using ``--node-data``) will automatically be included. These are
   often generated from the Augur commands ``traits``, ``clades`` or
   ``seqtraits``.
3. Traits present in the metadata file can be included by specifying
   them with ``metadata-color-by`` (e.g: ``--metadata-color-by country
   age host``). *(These must match column names of your metadata
   file.)*

The changes hopefully make things a little easier to use --
previously, if you had run ``augur clades``, you had to remember to add
``clade_membership`` to the config file, and if you'd run ``augur
seqtraits`` you had to add every resulting option. Now, they'll be
automatically included. If you don't want them as a coloring option,
don't pass in the files.


.. note:: You can't specify the title or type of a colouring option
          using just command-line - but ``export v2`` will make its
          best guess using the following rules: Excluding missing
          data, if a trait contains only 'True', 'False', 'Yes', 'No',
          '0' or '1', it will be set to 'boolean.' If it contains only
          numbers (integers and/or decimals), it will be set to
          'continuous.' Otherwise, it will be set as 'discrete.' If
          you want to have more control over how your trait is
          interpreted, you should use a config file (see :ref:`below
          <using-a-config-file-to-customise-the-visualisation>`).

Geographic Traits
*****************

Specify these traits using ``--geo-resolutions``, e.g.
``--geo-resolutions country region``. Previously these were defined by
the "geo" field in your v1 config file.

What's not possible to set without a config file
************************************************

The command line arguments cover everything you need to get a basic
run working in ``augur`` and ``auspice``. However, there are still
some features that offer more options or are only available when you
use a config file.

Currently, using command line arguments:

* It is not possible to set the default view options using only
  command-line arguments in ``export v2``. You can read more about the
  defaults (and how to change them using a config file)
  :ref:`here <display-defaults>`.
* When using ``export v2`` with only command-line arguments, every trait
  that's a coloring option and is either categorical or boolean will
  automatically be available to filter by. Find out how to specify
  what is a filter using a config file :ref:`here <filters>`.

-----------------------------------------------------------------

.. _using-a-config-file-to-customise-the-visualisation:

Using a config file to customise the visualisation
--------------------------------------------------

Traditionally you had to use an "Auspice config file" to customise the
visualisation. This is still available as an option, but you can now
choose between exporting using just the command-line, or using a
combination of the command-line and config file.

.. note:: Anything you can specify using the command-line arguments
          above can be done using a config file instead.

This section will detail the config file provided to ``augur export v2``
by the ``--auspice-config`` argument. The format of the new config file
**differs slightly** from previous versions of Augur. If you try to
use a previous version of the config file it *should mostly* still
work, but will print out warnings where keys have changed.


Config file priority
********************

It is important to remember that if you set an option both in the
config file *and* in the command line, **the command line option will
override the config file option**. For example, if you set ``"title"``
in your config file as "A Title About Apples", and then import this
config file using ``--auspice-config`` *and* use ``--title "Better Title
Befitting Bears"``, the title displayed by Auspice will be "Better
Title Befitting Bears". To use the one in the config file, don't use
``--title`` in the command line.

There are a couple of exceptions to this:

* There is no way to set default display views using command line
  only, so using "display_defaults" in your config file will set this.
* There is no way to modify the default filters displayed when using
  command line only, so using "filters" in your config file will set
  this.
* If you set color-by options in command-line using
  ``--metadata-color-by`` *and* pass in a config file, only the things
  listed in ``--metadata-color-by`` will be coloring options, but if
  they have a 'title' and 'type' set in the config file, these will be
  used.

Config file format
******************

The config file is a JSON file, and as such it's important that
everything in your config file is enclosed in one pair of curly
brackets. These can be on a separate line at the very top and very
bottom of your file. Syntax is important - if you are getting errors,
ensure all your brackets and quotation marks match up, and that commas
separate items in the same pair of brackets.

Export v2 config files are generally very similar to export v1, *but
there are a few changes*. They are explained in detail below, or you
can see :ref:`an example of converting a v1 config to v2 <updating-your-config-file>`.
For more details, see `the complete JSON schema for v2 config files
<https://github.com/nextstrain/augur/blob/-/augur/data/schema-auspice-config-v2.json>`__.

Here are the top-level keys of the config JSON in plain English:

title
,,,,,

The title to be displayed by Auspice, unchanged from previous versions
of the config file. E.g. ``"title": "Phylodynamics of my Pathogen"``.

maintainers
,,,,,,,,,,,

You can now have more than one maintainer associated with your run!
Specify one or as many maintainers as you wish via the following
structure (The ``url`` keys are optional):

.. code-block:: json

   "maintainers": [
     {"name": "Jane Doe", "url": "www.janedoe.com"},
     {"name": "Ravi Kupra", "url": "www.ravikupra.co.uk"}
   ]

Previously this was the "maintainer" field in your v1 config file and
used a different structure.

build-url
,,,,,,,,,

The build / repository URL to be displayed by Auspice, a new
functionality in ``augur export v2``, e.g. ``"build_url":
"https://github.com/nextstrain/zika"``. This is an optional field.

panels
,,,,,,

Optional and unchanged from previous versions of the config file, this
defines the panels that Auspice will display. If not set, Auspice will
by default try to show the tree, map, and entropy panels, if data is
available. Options are "tree", "map", "entropy", and "frequencies"
(e.g: ``"panels": ["tree", "map"]``).

.. note:: If you want to display the frequencies panel, you must
          specify "frequencies" *and* ensure a tip frequency file is
          available for ``auspice`` to access.

.. _colorings:

colorings
,,,,,,,,,

These are a list of the traits which Auspice should display as options
to color the tree & map. In previous versions of the config file this
was "color_options" and the current structure is similar, but
hopefully easier to understand!

For each trait you include, you can define:

* A required "key", which is used to lookup the values via node-data
  JSONs or other provided metadata.
* An optional "title" which will shown by Auspice when referring to
  this trait -- for instance you may have a trait called "ab1" which
  you want to show as "Age bracket 1" in the drop-down menus, legend,
  and filter.
* An optional, but highly recommended "type" which can be either
  'ordinal', 'boolean', 'continuous', or 'categorical'. If you don't
  provide a type, augur will try to guess it (see how it guesses
  :ref:`here <cli-traits>`).

Unless you want to change the name displayed, you *no longer* need to
include ``gt``, ``num_date``, ``clade_membership``, or ``augur
seqtraits`` output (like clade or drug resistance information) in your
config file - if that information is present, it will automatically be
included. To exclude it, don't pass in the corresponding file to
``--node-data``.

.. note:: Remember that if you are using ``--metadata-color-by`` on the
          command-line, only the traits given there will be color-by
          options! To include everything in your config file, don't
          use ``--metadata-color-by``, but include all traits you want
          as coloring options in "colorings" in the config file.

          Put another way, if a trait is listed in
          ``--metadata-color-by`` and not in the config, it will be
          included. If a trait is in the config but not in
          ``--metadata-color-by`` it will be excluded. If a trait is
          in both, but has ``"title"`` and ``"type"`` information in
          the config file, this information *will* be used by export v2.

In short, if using a config file and the command line, ensure
everything you want as a coloring option is in ``--metadata-color-by``.
You only need to also include it ``"colorings"`` in the config file if
you want to set the ``"title"`` and/or ``"type"``.

geo_resolutions
,,,,,,,,,,,,,,,

This specifies the geographical traits you want Auspice to use. You
can pass this in the same way as in the v1 config file, or you can now
specify a title to be displayed by option, using a slightly different
structure.

For example, for many users, these might be "country" and "region",
i.e. ``"geo_resolutions: ["country", "region"]``. If you want to give
them new titles, use the format ``"geo_resolutions": [{"key":
"country", "title": "Areas"}, {"key": "region", "title": "Global"}]``.

You can also mix the two, if you just want a title for one location:
``"geo_resolutions": [ {"key": "country", "title": "Areas"}, "region"]``

.. _filters:

filters
,,,,,,,

This specifies which traits you can filter by in Auspice. E.g.
``"filters": ["country", "region", "symptom", "age"]``. If you don't
include this option in your config file, all non-continuous traits
that are coloring options will be included as filters. If you don't
want any filter options, include this option with an empty list, i.e.
``"filters": []``. This is the same as the "filters" field in previous
config files, but the behavior has changed slightly.

.. _display-defaults:

display_defaults
,,,,,,,,,,,,,,,,

This allows you to specify the default view that users will see when
they visualise the data in Auspice. There are five options you can set
here -- note they are similar to those in the previous config files
but we have now standardised them to snake_case:

* ``geo_resolution`` - Sets which of the "geo_resolutions" should be
  shown. *Default: country*
* ``color_by`` - Sets what trait should be used for coloring. *Default: country*
* ``distance_measure`` - Sets whether tree branch lengths are in ``time``
  or ``divergence``. Options are ``num_date`` (time, default if available)
  or ``div`` (divergence).
* ``layout`` - Sets how the tree is visualised. Options are ``rect``
  (rectangular, default), ``radial``, ``unrooted``, and ``clock``,
  corresponding to the four options normally shown on the left in
  Auspice.
* ``map_triplicate`` - Sets whether the map is extended / wrapped
  around, which can be useful if transmissions are worldwide. Set to
  ``true`` or ``false``. *Default: false*

Config file examples
********************

Here is an example of how all of the above options would fit into a config file:

.. code-block:: json
   :linenos:

   {
     "title": "Phylodynamics of my Pathogen",
     "maintainers": [
       {"name": "Jane Doe", "url": "www.janedoe.com"},
       {"name": "Ravi Kupra", "url": "www.ravikupra.co.uk"}
     ],
     "build_url": "https://github.com/nextstrain/zika",
     "colorings": [
       {
         "key": "age",
         "title": "Host age",
         "type": "continuous"
       },
       {
         "key": "hospitalized",
         "type": "boolean"
       },
       {
         "key": "country",
         "type": "categorical"
       },
       {
         "key": "region",
         "type": "categorical"
       }
     ],
     "geo_resolutions": [
       {"key":"country", "title": "Areas"},
       "region"
     ],
     "panels": ["tree", "map"],
     "filters": ["country","region","symptom","age"],
     "display_defaults": {
       "color_by": "symptom",
       "geo_resolution": "region",
       "distance_measure": "div",
       "map_triplicate": "true"
     }
   }

If you want some examples of the new config files used in practice,
you can see some in these builds:

.. raw:: html

     <ul><li><span style="color:red"><em>TO DO</em></span></li></ul>



.. _updating-your-config-file:

Updating your config file
*************************

It's fairly easy to convert old export v1 config files to work with
export v2.

Here's an export v1 config file on the left, and an export v2 config
file on the right. We've tried to line them up to highlight the
differences:

+-----------------------------------------------+-----------------------------------------------------------------+
| Export v1 config:                             | Export v2 config:                                               |
|                                               |                                                                 |
| .. code-block:: json                          | .. code-block:: json                                            |
|    :linenos:                                  |    :linenos:                                                    |
|                                               |                                                                 |
|    {                                          |    {                                                            |
|      "title": "Phylodynamics of Virus A",     |      "title": "Phylodynamics of Virus A",                       |
|      "color_options": {                       |      "colorings": [                                             |
|        "gt": {                                |                                                                 |
|          "menuItem": "genotype",              |                                                                 |
|          "legendTitle": "Genotype",           |                                                                 |
|          "type": "discrete",                  |                                                                 |
|          "key": "genotype"                    |                                                                 |
|        },                                     |                                                                 |
|        "num_date": {                          |                                                                 |
|          "menuItem": "date",                  |                                                                 |
|          "legendTitle": "Sampling date",      |                                                                 |
|          "type": "continuous",                |                                                                 |
|          "key": "num_date"                    |                                                                 |
|        },                                     |                                                                 |
|        "age": {                               |        {                                                        |
|          "menuItem": "Host age",              |          "key": "age",                                          |
|          "legendTitle": "Host Age (years)",   |          "title": "Host age",                                   |
|          "type": "continuous",                |          "type": "continuous"                                   |
|          "key": "age"                         |        },                                                       |
|        },                                     |        {                                                        |
|        "host": {                              |          "key": "host",                                         |
|          "menuItem": "Animal",                |          "title": "Animal",                                     |
|          "legendTitle": "Animal",             |          "type": "categorical"                                  |
|          "type": "discrete",                  |        },                                                       |
|          "key": "host"                        |        {                                                        |
|        },                                     |          "key": "country",                                      |
|        "country": {                           |          "type": "categorical"                                  |
|          "type": "discrete"                   |        },                                                       |
|        },                                     |        {                                                        |
|        "region": {                            |          "key": "region",                                       |
|          "type": "discrete"                   |          "type": "categorical"                                  |
|        },                                     |        }                                                        |
|        "clade_membership":{                   |      ],                                                         |
|          "menuItem": "Clade",                 |                                                                 |
|          "legendTitle": "Clade",              |                                                                 |
|          "type": "discrete",                  |                                                                 |
|          "key": "clade_membership"            |                                                                 |
|        }                                      |                                                                 |
|      },                                       |                                                                 |
|      "geo": [                                 |      "geo_resolutions": [                                       |
|        "country",                             |          "country",                                             |
|        "region"                               |          "region"                                               |
|      ],                                       |        ],                                                       |
|      "maintainer": [                          |      "maintainers": [                                           |
|        "Hanna Kukk, Mohammad Fahir",          |        {"name": "Hanna Kukk", "url": "http://vamuzlab.org"},    |
|        "http://vamuzlab.org"                  |        {"name": "Mohammad Fahir", "url": "http://mfahir.co.uk"} |
|      ],                                       |      ],                                                         |
|      "filters": [                             |      "filters": [                                               |
|        "country", "region"                    |        "country", "region"                                      |
|      ],                                       |      ],                                                         |
|      "defaults": {                            |      "display_defaults": {                                      |
|        "layout": "unrooted",                  |        "layout": "unrooted",                                    |
|        "colorBy": "age"                       |        "color_by": "age"                                        |
|      }                                        |      }                                                          |
|    }                                          |    }                                                            |
+-----------------------------------------------+-----------------------------------------------------------------+

Vaccine choices
***************

In previous versions of augur, certain strains could be defined in the
config file as ``vaccine_choices`` (auspice would display this as a
cross over the tip in the tree). This functionality is now specified
via a node-data JSON (see the v6 release notes).

.. _prettifying-metadata-fields:

Prettifying metadata fields
---------------------------

In Auspice v1, we automatically 'prettified' many metadata values. For
example, a country value of 'new_zealand' would display as 'New
Zealand', and a metadata column called 'age_range' would display as
'Age Range'.

This worked well most of the time, but meant that users couldn't
intentionally keep underscores or lower-case values. It also meant we
had to detect exception cases like turning 'usa' into 'USA' rather
than 'Usa'.

In Auspice v2, all values are now displayed exactly as they arrive,
allowing users to ensure every gene and abbreviation displays just as
it should. However, this means that you should ensure your data looks
exactly how you'd like it to display - change any 'new_zealand's in
your metadata to 'New Zealand'!

Don't forget to also change them in any custom lat-long and/or
coloring files you are using. We've also become stricter about the
format of the files that pass in color and lat-long information.
Previously, it didn't matter if columns were separated by spaces or
tabs - now, they must be separated by tabs.

You can find out more about how to add
`custom coloring <https://docs.nextstrain.org/en/latest/guides/bioinformatics/colors.html>`__
and `lat-long <https://docs.nextstrain.org/en/latest/guides/bioinformatics/lat_longs.html>`__
values.

If you use the command ``parse`` to generate a metadata table from
fields in a fasta header, you can use the flag ``--prettify-fields`` to
apply some prettifying operations to specific metadata entries,
:doc:`see the documentation </usage/cli/parse>` for ``parse``.
