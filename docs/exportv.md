# How to Change 'Export' Versions

<span style="color:blue">

TO DO:
* **Update what's a colorby in CL to whatever we decide.**
* Update filter behaviour and explanation
* Explanation of config files and how/when overrides happen with CL. _Will be easier to write when we have resolved all questions ourselves_ (needs to include default view info)
* Decide how best to explain `--title` and `--maintainers` so that works both for CL and with snakemake _Current explanataion ok?_
* Confirm all options for `display_defaults`
* Do we need `--tree-name` (seems not always?)? In what circumstances? (Modify argument help description to reflect this)
* Do we support any of: `--output-sequence` `--reference` `--reference-translations` yet? In export only? In Auspice? If not, can we comment these out for the time being?

</span>

## What is `augur export v2`?

The `augur export` function is getting an upgrade! This is tied to an upgrade to Auspice (the visualization side of Nextstrain) to enable more functionality and flexibility. Your old 'version 1' (v1) JSON files will still work with the new version of `auspice`, but we recommend you switch to using 'version 2' (v2) JSONs from `export v2`, as this will become the new default. _Future versions of `auspice` may not support v1 JSON files._

### Why did we make this change?
* **Compactness**: Tree and Meta JSON files are now combined, so you only have to worry about one output file
* **Flexibility**: The new v2 JSONs allow us flexiblity to include more features and data, and will let us move towards getting in line with existing conventions like GFF and BibTex
* **Ease of use**: Users commonly got confused by the 'config' file. For basic runs you can now specify everything you need to see your data right in the command-line - no 'config' file needed! For more advanced exports, you can still specify a config file with more detail. (See [bottom of the page](#how-do-i-use-a-config-file-in-v2))

## **I just need my old run to work _right now_!**
*When you upgrade to the latest version of `augur`, `augur export` will no longer work.* But we understand you may not have time to make the change right this second, and that there's nothing more frustrating than having a run break right before a presentation or deadline!

If you want to keep using the old version _for now_, use `augur export v1` - everything else remains the same. 

To use the new version, use `augur export v2`. You'll need to make a few changes, but they are pretty simple, and you'll be future-proofing your runs! _(Future you will thank past you!)_

<br>


## Great! What do I need to do to move to `v2`?
You can always get a full overview of the arguments for export v2 with `augur export v2 --help`.

You can now choose between exporting using just the command-line or using a combination of the command-line and config file. Remember that **any command line options you use will override anything set in your config file**, if you are using both. 

We'll first cover what's the same/almost the same as in `export v1`, which are all things that must be passed in via command-line. Then we'll cover how to use [command-line options](#command-line-options), and finally how to use a [config file](#how-do-i-use-a-config-file-in-v2).


### What's the same:

You still pass in your tree, metadata, and node-data files with `--tree`, `--metadata`, and `--node-data` - just like in `export v1`.

Also just like in export v1, you can pass in files containing colors and latitute and longitude data using `--colors` and `--lat-longs`, respectively.

If you want to use a config file, you can pass this in with `--auspice-config`, just like in export v1.

### What's almost the same:

Instead of specifying two output files (`--output-tree` and `--output-meta`) you now just need to specify one with `--output-main`. 

For example, if your old files were `auspice/virus_AB_tree.json` and `auspice/virus_AB_meta.json`, you might want to call the single output `auspice/virus_AB.json` - or if you want to tell it apart from your v1 export, you might call it `auspice/virus_ABv2.json`. 

The URL for these, respectively, would be `...virus/AB` or `...virus/ABv2`.

In your metadata file, any column called `url` will be considered a link to that unique sequence in an online database (like Genbank), and will be the link attached to the Accession number. Any column called `paper_url` will be taken as a link associated to that author/title/journal. 

## Command-line Options

 #### General Display

* Specify the title of your run using `--title`. If running directly from the command line, put your title in quotes (ex: `--title "Phylodynamics of my Pathogen"`). If you are using snakemake and passing the value using `params`, you'll need to double-quote the title using single and double quotes. For example:
  ```
  params:
    title = "'Phylodynamics of my Pathogen'"
  shell:
    "augur export v2 --title {params.title} ..."
  ```

  _(Previously the "title" field in your config file.)_
<br>

* You can now have more than one maintainer associated with your run! Specify the maintainers with `--maintainers`. If running directly from the command line, put each maintainer in quotes (ex: `--maintainers "Jane Doe" "Ravi Kupra"`). If you are using snakemake and passing the value using `params`, you'll need to put the whole list in double quotes, and each person in single quotes. For example:
  ```
  params:
    maints = "'Jane Doe' 'Ravi Kupra'"
  shell:
    "augur export v2 --maintainers {params.maints} ..."
  ```
  You will need to use quotes in the same way even if you only have one maintainer!

   _(Previously the first part of the "maintainer" field in your config file.)_
<br>

* Specify the websites of maintainers with `--maintainer-urls`. Put them in the same order as the `--maintainers` so they link up with the right people! (ex: `--maintainer-urls www.janedoe.com www.ravikupra.co.uk`) **TO DO**

  _(Previously the second part of the "maintainer" field in your config file.)_
<br>

* If you want to specify what panels are visible, use `--panels`. By default, if the data is available, Auspice will show the tree, map, and entropy panels. 
You can specify "tree", "map", "entropy", and "frequencies". (ex: `--panels tree map entropy`) You must specify "frequencies" here _and_ supply a tip frequency file to `auspice` to display tip frequencies.

  _(Previously the "panels" field in the your config file.)_
<br>

* Specify how you'd like to position samples on the map using `--geography-traits`. For many users, these might be "country" and "region". (ex: `--geography-traits country region`)

  _(Previously the "geo" field in your config file.)_
<br>

  #### Traits

* Any traits that you have run with `augur traits` are now automatically included by `export v2` and available to color by! (These are often in a file called `traits_`(something)`.json`.) **TO DO**

  If you'd like to exclude any traits, you'll need to re-run `augur traits` without that trait. You can also do this with a config file (see [bottom of the page](#how-do-i-use-a-config-file-in-v2)).

  `export v2` will also automatically include date, genotype, and author information (if available) - you don't have to specify them!
<br>

* If you have extra meta-data that you'd like to include to color by, but
haven't run in `augur traits`, you can include it with `--extra-traits`.
Ensure you use the name of the column in the metadata file with the same
spelling and capitalization! 
(ex: `--extra-traits age host`)

  *(Previously, included traits were those things listed under "color_options" in your old config. `gt`, `num_date`, and `authors` don't need to be specified anymore - they'll be automatically included if present.)*
<br>

* If you don't provide a config file, `export v2` will try to 'guess' the type of the traits you include. Excluding missing data, if a trait contains only 'True', 'False', 'Yes', 'No', '0' or '1', it will be set to 'boolean.' If it contains only numbers (integers and/or decimals), it will be set to 'continuous.' Otherwise, it will be set as 'discrete.' If you want to have more control over how your trait is interpreted, you should use a config file (see [bottom of the page](#how-do-i-use-a-config-file-in-v2)).



### What's not possible in command-line only:
#### Default view
It is not possible to set the default view options using only command-line arguments in `export v2`. If the data is available, `auspice` will display your tree in rectangle view with branch lengths in time, colored by country, and plotted on the map by country. If time data is not available, it will set the branch lengths by divergence. If country data is not available as a geography trait, it will cycle through other options you passed in via `--geography-traits`. If country data is not available as a coloring option, it will cycle through other available colouring options. _What about `layout` option?_

If you would like to have more control over these options, you will need to include a config file (see [bottom of the page](#how-do-i-use-a-config-file-in-v2)). **TO DO** _Direct this link to a better place_

  _(Previously the "defaults" field in your config file.)_
<br>

#### Filter
When using `export v2` with only command-line arguments, every trait (both those from `augur traits` and passed in with `--extra-traits`), geographical trait, and the authors will automatically be available to filter by. Without using a config file, you can't exclude any of these from being filter options. **TO DO**


<span style="color:red">

_**TO DO**_

</span>


## How do I use a config file in `v2`?

### Why use a config file?
We have tried to make the command line options cover everything you need to get a run working in `augur` and `auspice`. However, there are still some features that offer more options or are only available when you use a config file. 

With a config file, you can:
* Specify exactly what you want as a filter option
* Set the default display view
* Give color-by traits more specific titles to display
* Specify color-by trait type
* _More advanced (flu?) options_ **TO DO**

### Config file priority
It is important to remember that if you set an option both in the config file _and_ in the command line, **the command line option will override the config file option**. For example, if you set `"title"` in your config file as "A Title About Apples", and then import this config file using `--auspice-config` _and_ use `--title "Better Title Befitting Bears"`, the title displayed by `auspice` will be "Better Title Befitting Bears". To use the one in the config file, simply don't use `--title` in the command line!

There are a couple of exceptions to this:
* There is no way to set default display views using command line only, so using `"display_defaults"` in your config file will set this
* There is no way to modify the default filters displayed when using command line only, so using `"filters"` in your config file will set this
* **???** _INCLUDE HOW THIS AFFECTS TRAITS DEPENDING ON HOW WE DECIDE THIS IN CL_
**TO DO**

### Using a Config File
You will always need to pass in some information, like your tree and metadata, via the command line. Jump back up to [What's the same](#what-s-the-same) and
[What's almost the same](#what-s-almost-the-same) to see what these are.

You'll also need to use `--auspice-config` to pass in the config file you'd like to use (this is the same as in export v1).

**Remember:** You don't have to include every field in your config file if you specify it in the command line (or are happy with the default) instead.

#### Config file format
Just like in export v1, the config file is a JSON-format file. In general, this means that it's divided into sections that then include more options and/or information in the curly (`{}`) or square (`[]`) brackets that follow.

**It is important that everything in your config file is enclosed in one pair of curly brackets.** These can be on a separate line at the very top and very bottom of your file. 

Indenting is not required for JSON files, but it can help make it easier to read your file. Syntax is important - if you are getting errors, check your brackets and quotation marks match up! Everything inside brackets is essentially a list, so needs commas between each item - but not after the last one!

_Export v2 config files are generally very simliar to export v1, but there are a few changes._ They are explained in detail below, or you can see [how to convert your v1 config to v2](#using-an-old-export-v1-config).

#### General Display
All of the options listed here go directly inside the main pair of curly brackets. See the [end of this section] for an example.

* Specify the title of your run using `"title"`. (ex: `"title": "Phylodynamics of my Pathogen"`)

  _(Same as the "title" field in your v1 config file.)_
<br>

* You can now have more than one maintainer associated with your run! Specify the runs maintainers and their websites using `"maintainers"` and listing the name and URL in pairs:
  ```
  "maintainers": [
    ["Jane Doe", "www.janedoe.com"], 
    ["Ravi Kupra","www.ravikupra.co.uk"]
  ]
  ```
  If you only have one maintiner, you still need to use the same format of two sets of square brackets: `"maintainers": [["Hanna Kukk", "www.hkukk.ee"]]`

  _(Previously the "maintainer" field in your v1 config file.)_
<br>

* To specify what panels are visible, use `"panels"`. By default, if the data is available, Auspice will show the tree, map, and entropy panels. 
You can specify "tree", "map", "entropy", and "frequencies". You must specify "frequencies" _and_ supply a tip frequency file to `auspice` to display tip frequencies.

  (ex: `"panels": ["tree", "map"]`) 

  _(Same as the "panels" field in the your v1 config file.)_
<br>

* Specify how you'd like to position samples on the map using `"geo"`. For many users, these might be "country" and "region". (ex: `"geo": [ "country", "region"]`)

  _(Same as the "geo" field in your v1 config file.)_
<br>

* Set what you would like to be able to filter by using `"filters"`. These must be traits present on your tree. (ex: `"filters": ["country", "region", "symptom", "age"]`)

  If you don't include this option in your config file, all non-continuous traits that are coloring options will be included as filters. If you don't want any filter options, include `"filters"` with no options (ex: `"filters": []`).

  _(Same as the "filters" field in your v1 config file.)_
<br>

* filters
* display options

* example


#### Traits

### Using an old (`export v1`) config
* Can this just go in as-is? Or what?

### Advanced use - second config file
* I don't think this exists yet?
 



