# How to Change 'Export' Versions

<span style="color:blue">

Emma's questions to be answered for this doc:
* Do we actually support config files now? _(Ensure this is up-to-date before we go live)_
* How do tip frequencies work in V2? Still separate file? Does 'panels' need to be specified to display this?
* Does panel specification work as I've described it? If you don't include a panel, is it definitely excluded? (Or what panels does this work for?)
* How do we specify default view, and does this work in `export v2`?
* When specifying `--title` in V2, why do we need to use quotes inside quotes? Can we do this better?
* How did (could?) we previously do multiple maintainers in old config file? (to include as comparison)

* Do we need `--tree-name` (seems not always?)? In what circumstances? (Modify argument help description to reflect this)
* Do we support any of: `--output-sequence` `--reference` `--reference-translations` yet? In export only? In Auspice? If not, can we comment these out for the time being?

</span>

## What is `augur export v2`?

The `augur export` function is getting an upgrade! This is tied to an upgrade to Auspice (the visualization side of Nextstrain) to enable more functionality and flexibility. Your old 'version 1' (v1) JSON files will still work with the new version of `auspice`, but we recommend you switch to using 'version 2' (v2) JSONs from `export v2`, as this will become the new default. _Future versions of `auspice` may not support v1 JSON files._

### Why did we make this change?
* **Compactness**: Tree and Meta JSON files are now combined, so you only have to worry about one output file
* **Flexibility**: The new v2 JSONs allow us flexiblity to include more features and data, and will let us move towards getting in line with existing conventions like GFF and BibTex
* **Ease of use**: Users commonly got confused by the 'config' file. For basic runs you can now specify everything you need to see your data right in the command-line - no 'config' file needed! For more advanced exports, you can still specify a config file with more detail. (See [bottom of the page](#how-do-i-use-a-config-file-in-v2)) _(coming soon)_

## **I just need my old run to work _right now_!**
*When you upgrade to the latest version of `augur`, `augur export` will no longer work.* But we understand you may not have time to make the change right this second, and that there's nothing more frustrating than having a run break right before a presentation or deadline!

If you want to keep using the old version _for now_, use `augur export v1` - everything else remains the same. 

To use the new version, use `augur export v2`. You'll need to make a few changes, but they are pretty simple, and you'll be future-proofing your runs! _(Future you will thank past you!)_

<br>


## Great! What do I need to do to move to `v2`?
You can always get a full overview of the arguments for export v2 with `augur export v2 --help`.

Here's how you can convert your v1 export and config file to  a **command-line-only** v2 export: 
_(See [bottom of page](#how-do-i-use-a-config-file-in-v2) for using a config file in v2)_


### What's the same:

You will still pass in your tree, metadata, and node-data files with `--tree`, `--metadata`, and `--node-data` - just like in `export v1`.

Also just like in export v1, you can pass in files containing colours and latitute and longitude data using `--colors` and `--lat-longs`, respectively.


### What's almost the same:

Instead of specifying two output files (`--output-tree` and `--output-meta`) you now just need to specify one with `--output-main`. 

For example, if your old files were `auspice/virus_AB_tree.json` and `auspice/virus_AB_meta.json`, you might want to call the single output `auspice/virus_AB.json` - or if you want to tell it apart from your v1 export, you might call it `auspice/virus_ABv2.json`. 

The URL for these, respectively, would be `...virus/AB` or `...virus/ABv2`.

### What's different:
* Specify the title of your run using `--title`. You'll need to use a bit of a strange quote system for this to work. For example: `--title '\'Phylodynamics of My Interesting Pathogen\''`  

  _(Previously the "title" field in your config file.)_
<br>

* You can now have more than one maintainer associated with your run! Specify the maintainers with `--maintainers`. Use double quotes for the  whole argument, and single-quotes to separate names (ex: `--maintainers "'Jane Doe' 'Ravi Kupra'"`)

   _(Previously the first part of the "maintainer" field in your config file.)_
<br>

* Specify the websites of maintainers with `--maintainer-urls`. Put them in the same order as the `--maintainers` so they link up with the right people! (ex: `--maintainer-urls 'www.janedoe.com www.ravikupra.co.uk'`)

  _(Previously the second part of the "maintainer" field in your config file.)_
<br>

* If you want to specify what panels are visible, use `--panels`. By default, if the data is available, Auspice will show the tree, map, and entropy panels. 
You can specify "tree", "map", "entropy", and "frequencies". _(???)_ (ex: `--panels "tree map entropy"`)

  _(Previously the "panels" field in the your config file.)_
<br>

* Specify how you'd like to locate samples on the map using `--geography-traits`. For many users, these might be "country" and "region". (ex: `--geography-traits "country region"`)

  _(Previously the "geo" field in your config file.)_
<br>

  #### Traits

* Any traits that you have run with `augur traits` are now automatically included by `export v2` and available to color by! (These are often in a file called `traits_`(something)`.json`.)

  If you'd like to exclude any traits, you'll need to re-run `augur traits` without that trait. You can also do this with a config file (see [bottom of the page](#how-do-i-use-a-config-file-in-v2)).

  `export v2` will also automatically include date, genotype, and author information (if available) - you don't have to specify them!
<br>

* If you have extra meta-data that you'd like to include to color by, but
haven't run in `augur traits`, you can include it with `--extra-traits`.
Ensure you use the name of the column in the metadata file with the same
spelling and capitalization! 
(ex: `--extra-traits "age host"`)

  *(Previously, included traits were those things listed under "color_options" in your old config. `gt`, `num_date`, and `authors` don't need to be specified anymore - they'll be automatically included if present.)*
<br>

* `export v2` will set traits as 'discrete' unless they contain only numbers, in which case they will be 'continuous.' If you want to have more control over how your trait is interpreted, you should use a config file (see [bottom of the page](#how-do-i-use-a-config-file-in-v2))



### What about filter?
When using `export v2` with only command-line arguments, every trait (both those from `augur traits` and passed in with `--extra-traits`), geographical trait, and the authors will automatically be available to filter by. Without using a config file, you can't exclude any of these from being filter options.


<span style="color:red">

_**TO DO**_

</span>


## How do I use a config file in `v2`?

### Why might I want to use a config file?
* exclude filter options
* make color-by traits pretty
* specify trait type (what do we support so far? discrete, continuous...)
* More advanced (flu?) options

### Config file priority
* Which overrides what and how, what must be included

### How one might start using a config
* Example of using a partial config to do simple things like traits and panels or filters

### Using an old (`export v1`) config
* Can this just go in as-is? Or what?

### Advanced use - second config file
* I don't think this exists yet?
 



