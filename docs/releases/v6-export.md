# Augur v6: exporting JSONs for auspice

<span style="color:orange">

TO DO:
* Check CL and config override stuff, expecially around colorings
* Decide how best to explain `--title` and `--maintainers` so that works both for CL and with snakemake _Current explanataion ok?_

</span>


While augur is a versitile bioinformatics toolkit, it was always developed in conjunction with [auspice](https://nextstrain.github.io/auspice), the interactive visualisation tool behind what you see on [nextstrain.org](https://nextstrain.org).
As such, most build pipelines finish with `augur export` which produces JSONs for auspice to visualise.

As auspice grows, the formats of the JSONs required changing (more on this in a bit) and, as such, augur's export function needed to change as well.
This page details how `augur export v2` works & what you need to do to start using it.

**Starting with augur v6, You must now specify which version of auspice you want to target when exporting from augur**.

---

* [Motivation behind changing JSON formats](#motivation-behind-changing-json-formats)
* [Compatability between augur & auspice versions](#compatability-between-augur-auspice-versions)
* [I just need my old run to work _right now_!](#i-just-need-my-old-run-to-work-right-now)
* [Terminology](#terminology)
* [Using command-line options to customise the visualisation](using-command-line-options-to-customise-the-visualisation)
* [Using a config file to customise the visualisation](#using-a-config-file-to-customise-the-visualisation)
---

## Motivation behind changing JSON formats

Augur v5 simply used the `augur export` command to produce JSONs (tree + meta JSON) for auspice v1. For this reason, well refer to these as "v1" JSONs.
With the release of auspice v2 [has now been released](https://nextstrain.github.io/auspice/releases/v2), we made some breaking changes to the JSON schema
This new JSON, which we refer to as "v2" (as it was designed as the input to auspice v2), was motivated by:

* **Compactness**: Tree and Meta JSON files are now combined, so you only have to worry about one output file
* **Flexibility**: The new v2 JSONs allow us flexiblity to include more features and data, and will let us move towards getting in line with existing conventions like GFF and BibTex
* **Ease of use**: Users commonly got confused by the 'config' file. For basic runs you can now specify everything you need to see your data right in the command-line - no 'config' file needed! For more advanced exports, you can still specify a config file with more detail. (See [bottom of the page](#config-file-options))

---
## Compatability between augur & auspice versions

Please see [this page](auspice-compatibility) for the most up-to-date compatibility between different augur & auspice versions.

> Please note that "v1" JSONs continue to work with auspice v2 -- we understand how important backwards compatability is.
We recommend switching to v2 JSONs as they have more features, are easier to work with, and one day future versions of auspice may no longer load v1 JSONs.


---
## I just need my old run to work _right now_!

We completely understand you may not have time to make the change right this second, and that there's nothing more frustrating than having a run break right before a presentation or deadline!

If you want to keep using the old version _for now_, replace `augur export` with `augur export v1` - everything else remains the same.
You can use auspice v1 or auspice v2 (see tables above).

To use the new version, use `augur export v2`.
You'll need to make a few changes, but they are pretty simple, and you'll be future-proofing your runs.
_(Future you will thank past you!)_

---
## What needs changing to use `augur export v2`?

> _P.S. You can always get a full overview of the arguments for export v2 with `augur export v2 --help`_

**What's the same**

You still pass in your tree, metadata, and node-data files with `--tree`, `--metadata`, and `--node-data` - just like in `export v1`.
Similarly, you can pass in files containing colors and latitute and longitude data using `--colors` and `--lat-longs`, respectively.
If you want to use a [config file](#config-file-options), you can pass this in with `--auspice-config`, but the format of this has changed (see below).

**Different outputs**

Instead of specifying two output files (`--output-tree` and `--output-meta`) you now just need to specify one with `--output`.
For example, if your old files were `auspice/virus_AB_tree.json` and `auspice/virus_AB_meta.json`, you might want to call the single output `auspice/virus_AB.json` - or if you want to tell it apart from your v1 export, you might call it `auspice/virus_ABv2.json`.


_TODO: the new option to export a reference sequence. Currently unused by Auspice_

**Command Line Options instead of (or in addition to) a config file**

The "auspice config" file defines a number of visualisation settings such as title, default displays and which colorings to use.
It has, however, always been a source of pain for a lot users.
You can now supply command line arguments instead of using a file if you wish (explained in more detail below).

---
## Terminology

#### Traits
Traits is the general term for certain data associated with nodes in the tree, for example "country" or "serotype" or "age".
These may have been inferred for internal nodes by the (confusingly named) `augur traits`, `augur clades` or others, or they may only be availiable for tips and provided by the metadata TSV file.

### Geographic Traits
Certain traits have a geographic interpretation, e.g. "country".
Auspice will attempt to display these traits on a map (and provide a drop-down to switch between them if there are more than one).

> _Make sure that these have corresponding entire in the lat-longs TSV file supplied to `export`._


---
## Using command-line options to customise the visualisation

As mentioned above, you can now replace most of the functionality of the "auspice config" file with command line options.
We hope that for most users this means the config file isn't necessary (but it's always there is you need its advanced functionality).

> Remember that generally _any command line options you use will override the same option in your config file_.


**Title**

Set the title displayed by auspice via `--title` (previously this was the "title" field in your v1 config file).
If running directly from the command line, put your title in quotes (ex: `--title "Phylodynamics of my Pathogen"`).
If you are using snakemake and passing the value using `params`, you'll need to double-quote the title using single and double quotes. For example:

```python
params:
  title = "'Phylodynamics of my Pathogen'"
shell:
  "augur export v2 --title {params.title} ..."
```


**Maintainers**

The maintainer(s) are displayed in the footer of auspice and may have associated links.
These can be specified with `--maintainers` and you can now have more than one maintainer associated with your run.
Previously this was set by the "maintainer" field in your v1 config file and was limited to a single entry.
If running directly from the command line, put each maintainer in quotes (ex: `--maintainers "Jane Doe" "Ravi Kupra"`). 
If you have a URL associated with a maintaner (completely optional), then you can add them like so:

```
--maintainers "Jane Doe <mailto:jane.doe@...>" "Ravi Kupra <https://github.com/ravikupra"
```

If you are using snakemake and passing the value using `params`, you'll need to put the whole list in double quotes, and each person in single quotes. For example:

```python
params:
  maints = "'Jane Doe' 'Ravi Kupra'"
shell:
  "augur export v2 --maintainers {params.maints} ..."
```
You will need to use quotes in the same way even if you only have one maintainer!


**Panels**

Auspice will, by default, try to show the tree, map, and entropy panels.
You can customise this with the `--panels` option, which was previously the "panels" field in the your v1 config file.
Options are "tree", "map", "entropy", and "frequencies". (e.g: `--panels tree map entropy`).

> If you want to display the frequencies panel, you must both specify "frequencies" here _and_ ensure a tip frequency file is available for `auspice` to access.


**Traits**

See above for the definition of traits.
Traits -- some automatically included, some defined on the command line -- will become coloring options in auspice.
The following rules are followed for which traits will be exported:

1. Genotype and date (if present) are always automatically included as coloring options - you don't need to include them.
(Previously these were "gt" and "numdate" in the "color_options" section of your v1 config file.)

2. Traits contained in the node-data JSONs handed to `augur export` will automatically be included.
These are often generated from the augur commands `traits`, `clades` or `seqtraits`.

3. Traits present in the metadata file can be included by specifying them with `metadata-color-by` (e.g: `--metadata-color-by country age host`).

The changes hopefully make things a little easier to use -- previously, if you had run `augur clades`, you had to remember to add `clade_membership` to the config file, and if you'd run `augur seqtraits` you had to add every resulting option.
Now, they'll be automatically included.
If you don't want them as a coloring option, just don't pass in the files!


> _Note: You can't specify the title or type of a colouring option using just command-line - but `export v2` will make its best guess using the following rules:_
_Excluding missing data, if a trait contains only 'True', 'False', 'Yes', 'No', '0' or '1', it will be set to 'boolean.'_
_If it contains only numbers (integers and/or decimals), it will be set to 'continuous.'_
_Otherwise, it will be set as 'discrete.'_
_If you want to have more control over how your trait is interpreted, you should use a config file (see [below](#config-file-options))._



**Geographic Traits**

Specify these traits using `--geo-resolutions`, e.g. `--geo-resolutions country region`.
Previously these were defined by the "geo" field in your v1 config file.


### What's not possible to set without a config file

We have tried to make the command line arguments cover everything you need to get a run working in `augur` and `auspice`.
However, there are still some features that offer more options or are only available when you use a config file.
Currently, using command line arguments:

* It is not possible to set the default view options using only command-line arguments in `export v2`.
* When using `export v2` with only command-line arguments, every trait that's a coloring option and is either categorical or boolean will automatically be available to filter by.



---
## Using a config file to customise the visualisation

Traditionally you had to use an "auspice config file" to customise the visualisation.
This is still available as an option but you can now choose between exporting using just the command-line, or using a combination of the command-line and config file.
Note that anything you can do with the command-line arguments detailed above can be done using a config file.

This section will detail the config file provided to `augur export v2` by the `--auspice-config` argument.
Note that the format of this file differs slightly to that used in previous versions of augur, mirroring the changes to the exported JSONs.
If you try to use a previous version of the config file it _should mostly_ still work, but will print out warnings where keys have changed etc.


#### Config file priority
It is important to remember that if you set an option both in the config file _and_ in the command line, _the command line option will override the config file option_.
For example, if you set `"title"` in your config file as "A Title About Apples", and then import this config file using `--auspice-config` _and_ use `--title "Better Title Befitting Bears"`, the title displayed by auspice will be "Better Title Befitting Bears".
To use the one in the config file, simply don't use `--title` in the command line!


There are a couple of exceptions to this:

* There is no way to set default display views using command line only, so using "display_defaults" in your config file will set this.
* There is no way to modify the default filters displayed when using command line only, so using "filters" in your config file will set this.
* If you set color-by options in command-line using `--metadata-color-by` _and_ pass in a config file, only the things listed in `--metadata-color-by` will be coloring options, but if they have a 'title' and 'type' set in the config file, these will be used.


#### Config file format

The config file is a JSON file, and as such it's important that everything in your config file is enclosed in one pair of curly brackets.
These can be on a separate line at the very top and very bottom of your file.
Syntax is important - if you are getting errors, ensure all your brackets and quotation marks match up, and that commas separate items in the same pair of brackets.


If you are familliar with JSON schemas, we have created one for this config file [here](https://github.com/nextstrain/augur/blob/v6/augur/data/schema-auspice-config-v2.json).


Here are the top-level keys of the config JSON in plain english:


**title**

The title to be displayed by Auspice, unchanged from previous versions of the config file.
E.g. `"title": "Phylodynamics of my Pathogen"`.

**maintainers**

You can now have more than one maintainer associated with your run!
Specify maintainers and their websites using `"maintainers"` and listing the name and URL in pairs:

```
"maintainers": [
  ["Jane Doe", "www.janedoe.com"],
  ["Ravi Kupra","www.ravikupra.co.uk"]
]
```

If you only have one maintainer, you still need to use the same format of two sets of square brackets: `"maintainers": [["Hanna Kukk", "www.hkukk.ee"]]`.
Previously this was the "maintainer" field in your v1 config file and used a different structure.

**panels**

Optional & unchanged from previous versions of the config file.
Defines the panels which auspice will display.
If not set then auspice will, by default, try to show the tree, map, and entropy panels.
Options are "tree", "map", "entropy", and "frequencies" (e.g: `"panels": ["tree", "map"]`).

> If you want to display the frequencies panel, you must both specify "frequencies" here _and_ ensure a tip frequency file is available for `auspice` to access.

**colorings**

Defined the traits which auspice should display as options to color the tree & map.
In previous versions of the config file this was "color_options" and the current structure is very similar but easier to understand!
For each trait defined here, you can define:
* an optional "title" which will shown by auspice when referring to this trait -- for instance you may have a trait called "ab1" which you want to show as "Age bracket 1". 
* an optional, but recommended "type" which must be one of  ordinal, boolean, continuous, or categorical. If you don't provide a type augur will try to guess it.


Unless you want to change the name displayed, you _no longer_ need to include `gt`, `num_date`, `clade_membership`, or `augur seqtraits` output (like clade or drug resistance information) in your config file - if that information is present, it will automatically be included. To exclude it, simply don't pass in the corresponding file to `--node-data`.

> _Remember that if you are using `--metadata-color-by` on the command-line, only the traits given there will be color-by options!_
_To include everything in your config file, don't use `--metadata-color-by`, rather include all traits you want as coloring options in "colorings" in the config file._
_Put another way, if a trait is listed in `--metadata-color-by` and not in the config, it will be included._
_If a trait is in the config but not in `--metadata-color-by` it will be excluded._
_If a trait is in both, but has `"title"` and `"type"` information in the config file, this information _will_ be used by export v2._
_In short, if using a config file and the command line, ensure everything you want as a coloring option is in `--metadata-color-by`._
_You only need to also include it `"colorings"` in the config file if you want to set the `"title"` and/or `"type"`._

**geo_resolutions**

Specify the geographical traits you want auspice to use.
For many users, these might be "country" and "region", i.e. `"geo_resolutions": [ {"key": country"}, {"key": "region"}]`.
This is _almost_ the same as the "geo" field in your v1 config file but the structure has changed slightly.


**filters**

Specify the traits you shich to show up as filters in auspice.
E.g. `"filters": ["country", "region", "symptom", "age"]`.
If you don't include this option in your config file, all non-continuous traits that are coloring options will be included as filters.
If you don't want any filter options, include this option with an empty list, i.e. `"filters": []`.
This is the same as the "filters" field in previous config files, but the behavior has changed slightly.


**display_defaults**

This allows you to specify the default view that users will see when they visualise the data in auspice.
There are five options you can set here -- note they are similar to those in the previous config files but we have now standardised them to snake_case:

* `geo_resolution` - Sets which of the "geo_resolutions" (see above) should be shown. 
* `color_by` - Sets what "coloring" to be used (see above).
* `distance_measure` - Sets whether tree branch lengths are in 'time' or 'divergence'.
Options are `num_date` (time, default if available) or `div` (divergence).
* `layout` - Sets how the tree is visualized.
Options are `rect` (rectangular, default), `radial`, `unrooted`, and `clock`, corresponding to the four options normally shown on the left in Auspice.
* `map_triplicate` - Sets whether the map is extended / wrapped around, which can be useful if transmissions are worldwide. Set to `true' or 'false'.



#### Config file examples

Here is an example of how all of the above options would fit into a config file _(Note this excludes trait color options, which is covered in the next section)_:
```json
{
  "title": "Phylodynamics of my Pathogen",
  "maintainers": [
    ["Jane Doe", "www.janedoe.com"],
    ["Ravi Kupra","www.ravikupra.co.uk"]
  ],
  "panels": ["tree", "map"],
  "colorings": {
    "age": {
      "title": "Host age",
      "type": "continuous"
    },
    "hospitalized": {
      "type": "boolean"
    },
    "country": {
      "type": "categorical"
    },
    "region": {
      "type": "categorical"
    }
  },
  "geo_resolutions": [
    {"key":"country"},
    {"key":"region"}
  ],
  "filters": ["country","region","symptom","age"],
  "display_defaults": {
    "color_by": "symptom",
    "geo_resolution": "region",
    "distance_measure": "div",
    "map_triplicate": "true"
  }
}
```


If you want some examples of the new config files used in practice, you can see some in these builds:
* _TO DO_

#### Updating your config file

It's fairly easy to convert old export v1 config files to work with export v2.

Here's an export v1 config file on the left, and an export v2 config file on the right. We've tried to line them up to highlight the differences:

<div style="overflow: hidden">
    <div id="column1" style="float:left; margin:0.5; width:50%;" markdown="1">
Export v1 config:
    <div class="highlight-default notranslate"><div class="highlight"><pre>
{
  "title": "Phylodynamics of Virus A",
  "color_options": {
    "gt": {
      "menuItem": "genotype",
      "legendTitle": "Genotype",
      "type": "discrete",
      "key": "genotype"
    },
    "num_date": {
      "menuItem": "date",
      "legendTitle": "Sampling date",
      "type": "continuous",
      "key": "num_date"
    },
    "age": {
      "menuItem": "Host age",
      "legendTitle": "Host Age (years)",
      "type": "continuous",
      "key": "age"
    },
    "host": {
      "menuItem": "Animal",
      "legendTitle": "Animal",
      "type": "discrete",
      "key": "host"
    },
    "country": {
      "type": "discrete"
    },
    "region": {
      "type": "discrete"
    },
    "clade_membership":{
      "menuItem": "Clade",
      "legendTitle": "Clade",
      "type": "discrete",
      "key": "clade_membership"
    }
  },
  "geo": [
    "country",
    "region"
  ],
  "maintainer": [
    "Hanna Kukk, Mohammad Fahir",
    "http://vamuzlab.org"
  ],
  "filters": [
    "country", "region"
  ],
  "defaults": {
    "layout": "unrooted",
    "colorBy": "age"
  }
}
    </pre></div></div>
    </div>
    <div id="column2" style="float:left; margin:0.5; width:50%;" markdown="1">
Export v2 config:
    <div class="highlight-default notranslate"><div class="highlight"><pre>
{
  "title": "Phylodynamics of Virus A",
  "colorings": {<br><br><br><br><br><br><br><br><br><br><br><br>
    "age": {
      "title": "Host age",<br>
      "type": "continuous"<br>
    },
    "host": {
      "title": "Animal",<br>
      "type": "categorical"<br>
    },
    "country": {
      "type": "categorical"
    },
    "region": {
      "type": "categorical"
    }<br><br><br><br><br><br>
  },
  "geo_resolutions": [
      {"key": "country"},
      {"key": "region"}
    ],
  "maintainers": [
    ["Hanna Kukk", "http://vamuzlab.org"],
    ["Mohammad Fahir", "http://mfahir.co.uk"]
  ],
  "filters": [
    "country", "region"
  ],
  "display_defaults": {
    "layout": "unrooted",
    "color_by": "age"
  }
}
    </pre></div></div>
    </div>
</div>