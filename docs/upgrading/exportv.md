# How to Change 'Export' Versions

<span style="color:blue">

TO DO:
* Check CL and config override stuff
* Decide how best to explain `--title` and `--maintainers` so that works both for CL and with snakemake _Current explanataion ok?_
* Are we including 'advanced config' options...?

</span>

## What is `augur export v2`?

The `augur export` function is getting an upgrade! This is tied to an upgrade to Auspice (the visualization side of Nextstrain) to enable more functionality and flexibility. Your old 'version 1' (v1) JSON files will still work with the new version of `auspice`, but we recommend you switch to using 'version 2' (v2) JSONs from `export v2`, as this will become the new default. _Future versions of `auspice` may not support v1 JSON files._

### Why did we make this change?
* **Compactness**: Tree and Meta JSON files are now combined, so you only have to worry about one output file
* **Flexibility**: The new v2 JSONs allow us flexiblity to include more features and data, and will let us move towards getting in line with existing conventions like GFF and BibTex
* **Ease of use**: Users commonly got confused by the 'config' file. For basic runs you can now specify everything you need to see your data right in the command-line - no 'config' file needed! For more advanced exports, you can still specify a config file with more detail. (See [bottom of the page](#config-file-options))

## **I just need my old run to work _right now_!**
*When you upgrade to the latest version of `augur`, `augur export` will no longer work.* But we understand you may not have time to make the change right this second, and that there's nothing more frustrating than having a run break right before a presentation or deadline!

If you want to keep using the old version _for now_, use `augur export v1` - everything else remains the same.

To use the new version, use `augur export v2`. You'll need to make a few changes, but they are pretty simple, and you'll be future-proofing your runs! _(Future you will thank past you!)_

<br>


## Great! What do I need to do to move to `v2`?
You can always get a full overview of the arguments for export v2 with `augur export v2 --help`.

You can now choose between exporting using just the command-line, or using a combination of the command-line and config file. Remember that generally **any command line options you use will override the same option in your config file**.

We'll first cover what's the same/almost the same as in `export v1`, which are all things that must be passed in via command-line. Then we'll cover how to use [command-line options](#command-line-options), and finally how to use a [config file](#config-file-options).


### What's the same:

You still pass in your tree, metadata, and node-data files with `--tree`, `--metadata`, and `--node-data` - just like in `export v1`.

Similarly, you can pass in files containing colors and latitute and longitude data using `--colors` and `--lat-longs`, respectively. If you want to use a [config file](#config-file-options), you can pass this in with `--auspice-config`.

### What's almost the same:

Instead of specifying two output files (`--output-tree` and `--output-meta`) you now just need to specify one with `--output-main`.

For example, if your old files were `auspice/virus_AB_tree.json` and `auspice/virus_AB_meta.json`, you might want to call the single output `auspice/virus_AB.json` - or if you want to tell it apart from your v1 export, you might call it `auspice/virus_ABv2.json`.

### How Coloring by Traits is Smarter
Previously config files always had to include `gt` and `num_date` to allow coloring by genotype and date, respectively. If you had run `augur clades`, you had to remember to add `clade_membership` to the config file, and if you'd run `augur seqtraits` you had to add every resulting option.

We've made things a little smarter! You never have to include `gt` or `num_date` in the command-line or config file - if the right data is sent to export v2, they'll automatically be included as a coloring option. The same is true for `augur clades` and `augur seqtraits` - if you pass the resulting JSON files in with `--node-data`, they'll be included as coloring options. _(If you don't want them as a coloring option, just don't pass in the files!)_

_Everything else_ you want to have as a coloring option must be specified either by `--metadata-color-by` or in the config file (what's specified in `--metadata-color-by` will override what's in the config file).

## Command-line Options
Remember that generally **any command line options you use will override the same option in your config file**.

 #### General Display

##### Title
* `--title` sets the title of your run
* If running directly from the command line, put your title in quotes (ex: `--title "Phylodynamics of my Pathogen"`). If you are using snakemake and passing the value using `params`, you'll need to double-quote the title using single and double quotes. For example:
  ```
  params:
    title = "'Phylodynamics of my Pathogen'"
  shell:
    "augur export v2 --title {params.title} ..."
  ```
* _(Previously the "title" field in your v1 config file.)_

##### Maintainers
* You can now have more than one maintainer associated with your run!
* Specify with `--maintainers`.
* If running directly from the command line, put each maintainer in quotes (ex: `--maintainers "Jane Doe" "Ravi Kupra"`). If you are using snakemake and passing the value using `params`, you'll need to put the whole list in double quotes, and each person in single quotes. For example:
  ```
  params:
    maints = "'Jane Doe' 'Ravi Kupra'"
  shell:
    "augur export v2 --maintainers {params.maints} ..."
  ```
  You will need to use quotes in the same way even if you only have one maintainer!
* _(Previously the first part of the "maintainer" field in your v1 config file.)_

##### Maintainer Websites
* Specify with `--maintainer-urls`.
* Put them in the same order as the `--maintainers` so they link up with the right people! (ex: `--maintainer-urls www.janedoe.com www.ravikupra.co.uk`)
* **TO DO** _We ok with this?_
* _(Previously the second part of the "maintainer" field in your v1 config file.)_

##### Panels
* `--panels` sets the visible panels.
* By default, if the data is available, Auspice will show the tree, map, and entropy panels.
* Options are "tree", "map", "entropy", and "frequencies". (ex: `--panels tree map entropy`) You must specify "frequencies" here _and_ supply a tip frequency file to `auspice` to display tip frequencies.
* _(Previously the "panels" field in the your v1 config file.)_

##### Geography
* Specify how you'd like to position samples on the map using `--geo-resolutions`.
* For many users, these might be "country" and "region". (ex: `--geo-resolutions country region`)
* _(Previously the "geo" field in your v1 config file.)_

#### Traits

* Genotype and date (if present) are always automatically included as coloring options - you don't need to include them.<br>  - _(Previously `gt` and `numdate` in your v1 config file)_
  <br><br>
* If you pass output from `augur clades` or `augur seqtraits` in using `--node-data`, that will also be automatically included as a coloring option <br>  - If you don't want them to be included, just don't pass the file in with `--node-data` <br>  - _(Previously things like `clade_membership`)_
<br><br>
* Pass in anything else from your metadata file you'd like to include with `--metadata-color-by` (ex: `--metadata-color-by country age host`) <br>  - _(Previously listed under "color_options" in your v1 config file)_
<br><br>
* You **can't** specify the title or type of a colouring option using just command-line - but `export v2` will make its best guess. <br>  - Excluding missing data, if a trait contains only 'True', 'False', 'Yes', 'No', '0' or '1', it will be set to 'boolean.' If it contains only numbers (integers and/or decimals), it will be set to 'continuous.' Otherwise, it will be set as 'discrete.' <br>  - If you want to have more control over how your trait is interpreted, you should use a config file (see [bottom of the page](#config-file-options)).




### What's not possible in command-line only:
#### Default view
* It is not possible to set the default view options using only command-line arguments in `export v2`.

* By default, `auspice` will display your tree in rectangle view, with branch lengths, color, and map position based on the available data. You can read more about the defaults [here](#id6).

* If you would like to have more control over these options, you will need to include a config file. You can find more information [here](#id6).

* _(Previously the "defaults" field in your v1 config file.)_
<br>

#### Filter
* When using `export v2` with only command-line arguments, every trait that's a coloring option and is either categorical or boolean will automatically be available to filter by.

* If you'd like to specify exactly what traits are listed as filters, you'll need to use a config file. You can find more information [here](#filters).

*  _(The same as the "filters" field in your v1 config file.)_



## Config File Options

### Why use a config file?
We have tried to make the command line arguments cover everything you need to get a run working in `augur` and `auspice`. However, there are still some features that offer more options or are only available when you use a config file.

With a config file, you can:
* Specify exactly what you want as a filter option
* Set the default display view
* Give color-by traits more specific titles
* Specify color-by trait type
* _More advanced (flu?) options_ **TO DO**

### Config file priority
It is important to remember that if you set an option both in the config file _and_ in the command line, **the command line option will override the config file option**. For example, if you set `"title"` in your config file as "A Title About Apples", and then import this config file using `--auspice-config` _and_ use `--title "Better Title Befitting Bears"`, the title displayed by `auspice` will be "Better Title Befitting Bears". To use the one in the config file, simply don't use `--title` in the command line!

There are a couple of exceptions to this:
* There is no way to set default display views using command line only, so using `"display_defaults"` in your config file will set this
* There is no way to modify the default filters displayed when using command line only, so using `"filters"` in your config file will set this
* If you set color-by options in command-line using `--metadata-color-by` _and_ pass in a config file, only the things listed in `--metadata-color-by` will be coloring options, but if they have a 'title' and 'type' set in the config file, these will be used.


### Using a Config File
* You will always need to pass in some information, like your tree and metadata, via the command line. Jump back up to [What's the same](#what-s-the-same) and
[What's almost the same](#what-s-almost-the-same) to see what these are.

* You'll also need to use `--auspice-config` to pass in the config file you'd like to use (this is the same as in export v1).

* **Remember:** You don't have to include every field in your config file - you can specify it in the command line (or use the default default) instead, if you prefer.

#### Config file format
Just like in export v1, the config file is a JSON-format file. _It is important that everything in your config file is enclosed in one pair of curly brackets._ These can be on a separate line at the very top and very bottom of your file.

Syntax is important - if you are getting errors, ensure all your brackets and quotation marks match up, and that commas separate items in the same pair of brackets.

_Export v2 config files are generally very simliar to export v1, but there are a few changes._ They are explained in detail below, or you can see [an example of converting a v1 config to v2](#using-an-old-export-v1-config).

#### General Display
All of the options listed here go directly inside the main pair of curly brackets. See the [end of this section](#example) for an example.

##### Title
* Specify the title of your run using `"title"`.
* Ex: `"title": "Phylodynamics of my Pathogen"`
* _(Same as the "title" field in your v1 config file.)_

##### Maintainers
* You can now have more than one maintainer associated with your run!
* Specify maintainers and their websites using `"maintainers"` and listing the name and URL in pairs:
  ```
  "maintainers": [
    ["Jane Doe", "www.janedoe.com"],
    ["Ravi Kupra","www.ravikupra.co.uk"]
  ]
  ```
* If you only have one maintainer, you still need to use the same format of two sets of square brackets: `"maintainers": [["Hanna Kukk", "www.hkukk.ee"]]`
* _(Previously the "maintainer" field in your v1 config file.)_

##### Panels
* Use `"panels"` to specify visible panels.
* By default, if the data is available, Auspice will show the tree, map, and entropy panels.
* Options are "tree", "map", "entropy", and "frequencies". You must specify "frequencies" _and_ supply a tip frequency file to `auspice` to display tip frequencies.
* Ex: `"panels": ["tree", "map"]`
* _(Same as the "panels" field in the your v1 config file.)_

##### Geography
* Specify how you'd like to position samples on the map using `"geo_resolutions"`.
* For many users, these might be "country" and "region". (ex: `"geo_resolutions": [ {"key": country"}, {"key": "region"}]`)
* _(Almost the same as the "geo" field in your v1 config file.)_

##### Filters
* Set what you would like to be able to filter by using `"filters"`.
* These must be traits present on your tree.
* Ex: `"filters": ["country", "region", "symptom", "age"]`
* If you don't include this option in your config file, all non-continuous traits that are coloring options will be included as filters. If you don't want any filter options, include `"filters"` with no options (ex: `"filters": []`).
* _(Same as the "filters" field in your v1 config file.)_

##### Default View
* You can specify the default view that users will see when they load your run by using `"display_defaults"`.
* If you do not change the options here, `auspice` will revert to the defaults listed below:

  There are five options you can set here:
  * `geo_resolution` - Sets which `geo` option is used to position data on the map. Default is `country` - if not available, other `geo` options will be tried.
  * `color_by` - Sets what coloring trait the tree should be colored by. Must be an available coloring trait. Default is `country` - if not available, other coloring options will be tried.
  * `distance_measure` - Sets whether tree branch lengths are in 'time' or 'divergence'. Default is `num_date` (time), if available. Options are `num_date` (time) or `div` (divergence).
  * `layout` - Sets how the tree is visualized. Default is `rect` (rectangle). Options are `rect`, `radial`, `unrooted`, and `clock`, corresponding to the four options normally shown on the left in Auspice.
  * `map_triplicate` - Sets whether the map is extended / wrapped around, which can be useful if transmissions are worldwide. Set to 'true' or 'false'.
  <br>
* _(Previously the "defaults" field in your config file. Note the spelling of the five options has changed slightly!)_

##### Example
Here is an example of how all of the above options would fit into a config file _(Note this excludes trait color options, which is covered in the next section)_:
  ```
  {
    "title": "Phylodynamics of my Pathogen",
    "maintainers": [
      ["Jane Doe", "www.janedoe.com"],
      ["Ravi Kupra","www.ravikupra.co.uk"]
    ],
    "panels": ["tree", "map"],
    "geo_resolutions": [ {"key":"country"}, {"key":"region"}],
    "filters": [
      "country", "region", "symptom", "age"
    ],
    "display_defaults": {
      "color_by": "symptom",
      "geo_resolution": "region",
      "distance_measure": "div",
      "map_triplicate": "true"
    }
  }
  ```

#### Traits
**Remember that if you are using `--metadata-color-by` on the command-line, only the traits given there will be color-by options! To include everything in your config file, don't use `--metadata-color-by`.

In the export v1 config file, `"color_options"` was used to specify and describe traits that you wanted as colouring options. This is now replaced by `"colorings"`, and works in a very similar way.

As in export v1, the `"colorings"` section of the config file will hold all the trait information within its curly brackets (see [the example](#example)). You should provide traits using the column name in the metadata.

If you wish, you can provide more information about the trait using `"title"` (what should display in the drop-down menu in Auspice - previously `"menuItem"` in export v1) and `"type"` (should the trait be treated as ordinal, boolean, continuous, or categorical). If you don't provide these, export v2 will simply use the trait name as provided and try to guess the type. The export v1 options `"legendTitle"` and `"key"` are no longer used.

Unless you want to change the name displayed, you _no longer_ need to include `gt`, `num_date`, `clade_membership`, or `augur seqtraits` output (like clade or drug resistance information) in your config file - if that information is present, it will automatically be included. To exclude it, simply don't pass in the corresponding file to `--node-data`.

##### Config file only
To specify coloring options using only a config file, **do not** use `--metadata-color-by` in the command line. Then include all traits you want as coloring options in `"colorings"` in the config file. As mentioned above, you can include `"title"` and `"type"` information, or not.

##### Config file and Command Line
**Using `--metadata-color-by` on the command line will override what's included in `"colorings"` in the config file to decide what will be a coloring option.** If a trait is listed in `--metadata-color-by` and not in the config, it will be included. If a trait is in the config but not in `--metadata-color-by` it will be excluded. _However_, if a trait is in both, but has `"title"` and `"type"` information in the config file, this information _will_ be used by export v2.

In short, if using a config file and the command line, ensure everything you want as a coloring option is in `--metadata-color-by`. You only need to also include it `"colorings"` in the config file if you want to set the `"title"` and/or `"type"`.

##### Example
Here's an example of how display and trait options (for `age`, `hospitalized`, `country`, and `region`) might fit together in a config file:
  ```
  {
    "title": "Phylodynamics of my Pathogen",
    "colorings": {
      "age": {
        "title": "Host age",
        "type": "continuous"
      },
      "hospitalized": {
        "type": "boolean"
      },
      "country": {
      },
      "region": {
      }
    },
    "maintainers": [
      ["Jane Doe", "www.janedoe.com"],
      ["Ravi Kupra","www.ravikupra.co.uk"]
    ],
    "geo_resolutions": [ {"key":"country"}, {"key":"region"}],
    "filters": [ "country", "region" ],
    "display_defaults": {
      "geo_resolution": "region",
    }
  }
  ```

### Using an old (`export v1`) config
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
    "country", "region"
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
      {"key":"country"},
      {"key":"region"}
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

#### What changed?

* `title` and `filters` are unchanged
* `"color_options"` is now `"colorings"`. In this section:
  * `gt`, `num_date`, and `clade_membership` no longer have to be included - they are automatically included if the necessary information is present
  * `legendTitle` and `key` are gone
  * `menuItem` is now `title`
  * The type `discrete` is now `categorical`
* `geo` is now `geo_resolutions` and uses a new format
* `maintainer` is now `maintainers` and uses the new format
* `defaults` is now `display_defaults` and within its options, `colorBy` is now `color_by`


### Advanced use - second config file
**TO DO**
* I don't think this exists yet?
_Are we doing this?_




