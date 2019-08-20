# Breaking Changes to `augur`

This is a list of known 'breaking changes' to `augur`. These are changes that require users to modify existing scripts or calls to `augur` functions because the old usage is no longer supported. 'Deprecated' changes, which work for the moment but will no longer be supported in future, are also included. 

_Note this list only includes changes to `augur` code and may not cover breaking changes introduced due to changes in dependancy packages. We try to support these too, so please [open an issue](https://github.com/nextstrain/augur/issues/new) if you think you've found something we didn't catch._

This list will **only** be helpful if you are _sure_ the command used to work and only stopped working after you upgraded your `augur` installation. If your `augur` installation is from earlier than 2019, or if this list doesn't solve your problem, run `--help` for the command you're using and use that information to try re-writing your call.

Click on the `augur` function that used to work, and we'll try to get you up and running ASAP!

* [ancestral](#ancestral)
* [export](#export)

## ancestral

### `--output` argument

* **Possible error/warning messages:**<br>
  > > "WARNING: the `--output` flag will be deprecated in the next major augur release. Use `--output-node-data` instead."

  > > augur: error: unrecognized arguments: `--output`

* **Solution:**<br>
  To specify the JSON file to write ancestral mutations/sequences to, use the argument `--output-node-data` instead of `--output`.<br><br>

* **Explanation:**<br>
  `ancestral` now supports outputting a Fasta file of reconstructed ancestral sequences (for Fasta-input runs - this was already supported in VCF-format for VCF-input runs). Users can ask for this output and specify a file name using `--output-sequences`. Since there are now two types of output from `ancestral`, the arguments have become more descriptive.<br><br>


## export

### Version error

* **Possible error/warning messages:**<br>
  > > augur export: error: the following arguments are required: Augur export now needs you to define the JSON version you want, e.g. `augur export v2`.

* **Solution:**<br>
  We've upgraded `augur export`. Find out how to adjust your `augur export` call to work with the latest version [using our handy guide](exportv.md).<br><br>
