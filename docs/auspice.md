# How Augur Export and Auspice Versions Connect

From mid-2019, we're changing the format of the output file produced by `augur export` to allow us more flexibility going forward. This lines up with a big update to `auspice`, which also opens up exciting new possibilities. 

We also recognise that we may need to adjust things again in the future. By setting up versions now, we'll be able to continue to expand `augur` and `auspice` while minimizing disruption to users.

You can get more information on how to move to using `export v2` in `augur` 6.0 in our [handy guide](exportv.md).

## Why Does Compatibility Matter?

To get great visualizations, the final step of the `augur` pipeline, `augur export` produces output as (a) JSON file(s)*. This/these file(s) are then read by `auspice` to produce trees, graphs, and maps. 

*_This compatibility only refers to the JSON file(s) produced by `augur export`, not those produced by other `augur` steps._

Different `augur` versions will produce different JSON output formats, which in turn will work with different `auspice` versions. We aim to support backwards compatibility as much as possible, and of course, if you're you're using the latest `augur`, JSON, and `auspice` versions, everything will work smoothly!

Whenever possible, we recommend updating your run to work with the latest `augur` and outputing the most recent JSON format to view with the latest `auspice` - it future-proofs your work and ensures you'll benefit from all the latest features and improvements.

But, we know that for a variety of reasons, you can't always do this straight away. To help keep your research flowing, here are some tables of how different versions of `augur`, export JSONs, and `auspice` are compatible.

## Compatibility Tables

### Augur

| Augur Version | Produces JSON | Works with Auspice  |
| ------------- |---------------| --------------------|
| 5.x           | v1 (via `export`)   | 1.x or 2.0          |
| 6.0           | v1 (via `export v1`) <br> v2 (via `export v2`) | 1.x or 2.0 |


### `export` JSON versions

| JSON Version | Produced by Augur Version | Works with Auspice  |
| -------------|---------------| --------------------|
| v1           | 5.x (via `export`) <br> 6.0 (via `export v1`)| 1.x or 2.0 |
| v2           | 6.0 (via `export v2`) | 2.0 |

### Auspice

| Auspice Version | Takes JSON | Produced by Augur Version  |
| -------------|---------------| --------------------|
| 1.x           | v1        | 5.x (via `export`) <br> 6.0 (via `export v1`) |
| 2.0           | v1 <br> v2 | 5.x (via `export`) <br> 6.0 (via `export v1`) <br> 6.0 (via `export v2`) |
