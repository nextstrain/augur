# CHANGELOG

## __NEXT__

### Features

* A new command, `augur merge`, now allows for generalized merging of two or more metadata tables. [#1563][] (@tsibley)
* Two new commands, `augur read-file` and `augur write-file`, now allow external programs to do i/o like Augur by piping from/to these new commands.  They provide handling of compression formats and newlines consistent with the rest of Augur. [#1562][] (@tsibley)
* A new debugging mode can be enabled by setting the `AUGUR_DEBUG` environment variable to `1` (or another truthy value).  Currently the only effect is to print more information about handled (i.e. anticipated) errors.  For example, stack traces and parent exceptions in an exception chain are normally omitted for handled errors, but setting this env var includes them.  Future debugging and troubleshooting features, like verbose operation logging, will likely also condition on this new debugging mode. [#1577][] (@tsibley)

### Bug Fixes

* Embedded newlines in quoted field values of metadata files read/written by many commands, annotation files read by `augur curate apply-record-annotations`, and index files written by `augur index` are now properly handled. [#1561][] [#1564][] (@tsibley)

[#1561]: https://github.com/nextstrain/augur/pull/1561
[#1562]: https://github.com/nextstrain/augur/pull/1562
[#1563]: https://github.com/nextstrain/augur/pull/1563
[#1564]: https://github.com/nextstrain/augur/pull/1564
[#1577]: https://github.com/nextstrain/augur/pull/1577



## 25.2.0 (24 July 2024)

### Features

* export v2: we now limit numerical precision on floats in the JSON. This should not change how a dataset is displayed / interpreted in Auspice but allows the gzipped & minimised JSON filesize to be reduced by around 30% (dataset-dependent). [#1512][] (@jameshadfield)
* traits, export v2: `augur traits` now reports all confidence values above 0.1% rather than limiting them to the top 4 results. There is no change in the eventual Auspice dataset as `augur export v2` will still only consider the top 4. [#1512][] (@jameshadfield)
* curate: Excel (`.xlsx` and `.xls`) and OpenOffice (`.ods`) spreadsheet files are now also supported as metadata inputs (`--metadata`).  The first sheet in the workbook is read as tabular data.  [#1550][] (@tsibley)

### Bug Fixes

* titers sub: Fixes a bug where antigenic weights were assigned to branches for substitutions in the incorrect order of `<derived allele><position><ancestral allele>` instead of `<ancestral allele><position><derived allele>`. [#1555][] (@huddlej)

[#1512]: https://github.com/nextstrain/augur/pull/1512
[#1550]: https://github.com/nextstrain/augur/pull/1550
[#1555]: https://github.com/nextstrain/augur/pull/1555


## 25.1.1 (15 July 2024)

### Bug Fixes

* curate parse-genbank-location: Fix a bug where a mix of empty and populated location-field values would result in inconsistent fields in the output NDJSON [#1531][](@genehack)

[#1531]: https://github.com/nextstrain/augur/pull/1531


## 25.1.0 (11 July 2024)

### Features

* Support xopen major version 2. Deprecate v1. Schedule for removal around November 2024. [#1532][] (@corneliusroemer)
* Support networkx major version 3. [#1534][] (@corneliusroemer)

[#1532]: https://github.com/nextstrain/augur/pull/1532
[#1534]: https://github.com/nextstrain/augur/pull/1534

## 25.0.0 (10 July 2024)

### Major changes

* curate format-dates: Raises an error if provided date field does not exist in records. [#1509][] (@joverlee521)
* All curate subcommands: Verifies all input records have the same fields and raises an error if a record does not have matching fields. [#1518][] (@joverlee521)

### Features

* Added a new sub-command `augur curate apply-geolocation-rules` to apply user curated geolocation rules to the geolocation fields in a metadata file. Previously, this was available as a script within the nextstrain/ingest repo. [#1491][] (@victorlin)
* Added a default color for the "Asia" region that will be used in `augur export` is no custom colors are provided. [#1490][] (@joverlee521)
* Added a new sub-command `augur curate apply-record-annotations` to apply user curated annotations to existing fields in a metadata file. Previously, this was available as a `merge-user-metadata` in the nextstrain/ingest repo. [#1495][] (@joverlee521)
* Added a new sub-command `augur curate abbreviate-authors` to abbreviate lists of authors to "<first author> et al." Previously, this was avaliable as the `transform-authors` script within the nextstrain/ingest repo. [#1483][] (@genehack)
* Added a new sub-command `augur curate parse-genbank-location` to parse the `geo_loc_name` field from GenBank reconds. Previously, this was available as the `translate-genbank-location` script within the nextstrain/ingest repo. [#1485][] (@genehack)
* curate format-dates: Added defaults to `--expected-date-formats` so that ISO 8601 dates (`%Y-%m-%d`) and its various masked forms (e.g. `%Y-XX-XX`) are automatically parsed by the command. [#1501][] (@joverlee521)
* Added a new sub-command `augur curate transform-strain-name` to filter strain names based on matching a regular expression. Previously, this was available as the `transform-strain-names` script within the nextstrain/ingest repo. [#1514][] (@genehack)
* Added a new sub-command `augur curate rename` to rename field / column names. Previously, a similar version was available as the `transform-field-names` script within the nextstrain/ingest repo however the behaviour is slightly changed here. [#1506][] (@jameshadfield)

### Bug Fixes

* filter: Improve speed of checking duplicates in metadata, especially for large files. [#1466][] (@victorlin)
* curate: Stop adding double quotes to the metadata TSV output when field values have internal quotes. [#1493][] (@joverlee521)
* curate format-dates: Mask empty date values as `XXXX-XX-XX` to represent unknown dates. [#1509][] (@joverlee521)

[#1466]: https://github.com/nextstrain/augur/pull/1466
[#1490]: https://github.com/nextstrain/augur/pull/1490
[#1491]: https://github.com/nextstrain/augur/pull/1491
[#1493]: https://github.com/nextstrain/augur/pull/1493
[#1495]: https://github.com/nextstrain/augur/pull/1495
[#1501]: https://github.com/nextstrain/augur/pull/1501
[#1509]: https://github.com/nextstrain/augur/pull/1509
[#1514]: https://github.com/nextstrain/augur/pull/1514
[#1518]: https://github.com/nextstrain/augur/pull/1518
[#1506]: https://github.com/nextstrain/augur/pull/1506

## 24.4.0 (15 May 2024)

### Features

* All commands: Allow repeating an option that takes multiple values. Previously, if multiple option flags were specified (e.g. `--exclude-where 'region=A' --exclude-where 'region=B'`), only the last one was used. Now, all values are used. [#1445][] (@victorlin)
* ancestral, translate: output node data files are now validated. The argument `--validation-mode` is added which controls this behaviour (default: error). This argument also controls validation of the input node-data file (ancestral only). [#1440][] (@jameshadfield)
* export: Updated default latitudes and longitudes for geography traits. This only applies if you are **not** using `--lat-longs` to override the built in mappings. [#1449][] (@trvrb)

### Bug Fixes

* validation: we no longer exit with a non-zero exit code when the requested validation mode is "warn" [#1440][] (@jameshadfield)
* validation: we no longer perform any validation when the requested validation mode is "skip" [#1440][] (@jameshadfield)
* filter: Send all log messages to `stderr`. This allows output to be written to `stdout` (e.g. `--output-strains /dev/stdout`). [#1459][] (@victorlin)

[#1440]: https://github.com/nextstrain/augur/pull/1440
[#1445]: https://github.com/nextstrain/augur/pull/1445
[#1449]: https://github.com/nextstrain/augur/pull/1449
[#1459]: https://github.com/nextstrain/augur/pull/1459

## 24.3.0 (18 March 2024)

### Features

* filter: Added a new option `--max-length` to filter out sequences that are longer than a certain amount of base pairs. [#1429][] (@victorlin)
* parse: Added support for environments that use pandas 2.x. [#1436][] (@emollier, @victorlin)

### Bug Fixes

* filter: Updated docs with an example of tiered subsampling. [#1425][] (@victorlin)
* export: Fixes bug [#1433] introduced in v23.1.0, that causes validation to fail when gene names start with `nuc`, e.g. `nucleocapsid`. [#1434][] (@corneliusroemer)
* import: Fixes bug introduced in v24.2.0 that prevented `import beast` from running. [#1439][] (@tomkinsc)
* translate, ancestral: Compound CDS are now exported as segmented CDS and are now viewable in Auspice.  [#1438][] (@jameshadfield)

[#1425]: https://github.com/nextstrain/augur/pull/1425
[#1429]: https://github.com/nextstrain/augur/pull/1429
[#1433]: https://github.com/nextstrain/augur/issues/1433
[#1434]: https://github.com/nextstrain/augur/pull/1434
[#1436]: https://github.com/nextstrain/augur/pull/1436
[#1438]: https://github.com/nextstrain/augur/pull/1438
[#1439]: https://github.com/nextstrain/augur/pull/1439

## 24.2.3 (23 February 2024)

### Bug Fixes

* filter: Updated the help and report text of `--min-length` to explicitly state that the minimum length filter only counts standard nucleotide characters A, C, G, or T (case-insensitive). This has been the behavior since version 3.0.3.dev1, but has never been explicitly documented. [#1422][] (@joverlee521)
* frequencies: Fixed a bug introduced in 24.2.0 and 24.1.0 that prevented `--regions` from working when providing regions other than the default "global" region. [#1424]

[#1422]: https://github.com/nextstrain/augur/pull/1422
[#1424]: https://github.com/nextstrain/augur/pull/1424

## 24.2.2 (16 February 2024)

### Bug Fixes

* filter: In versions 24.2.0 and 24.2.1, `--query` stopped working in cases where internal optimizations added in version 24.2.0 failed to parse the columns from the query. It now falls back to non-optimized behavior that allows queries to work. [#1418][] (@victorlin)
* filter: Handle backtick quoting in internal optimizations of `--query`. [#1417][] (@victorlin)

[#1417]: https://github.com/nextstrain/augur/pull/1417
[#1418]: https://github.com/nextstrain/augur/pull/1418

## 24.2.1 (14 February 2024)

### Bug Fixes

* frequencies: Fixed a bug introduced in 24.2.0 that prevented `--method diffusion` from working alongside `--tree`. [#1412][] (@victorlin)

[#1412]: https://github.com/nextstrain/augur/issues/1412

## 24.2.0 (12 February 2024)

### Features

* filter: Added a new option `--query-columns` that allows specifying what columns are used in `--query` along with the expected data types. If unspecified, automatic detection of columns and types is attempted. [#1294][] (@victorlin)
* `augur.io.read_metadata`: A new optional `columns` argument allows specifying a subset of columns to load. The default behavior still loads all columns, so this is not a breaking change. [#1294][] (@victorlin)
* `augur parse`: A new optional `--output-id-field` argument allows the user to select any ID field for the produced FASTA file (e.g. 'accession' instead of 'name' or 'strain'). [#1403][] (@j23414)
  * When no `--output-id-field` is given and the data has both `name` and `strain` fields, continue to preferentially use `name` over `strain` as the sequence ID field; but, throw a deprecation warning that the order will be switched to prefer `strain` over `name` in the future to be consistent with the rest of Augur.
  * Added entry to [DEPRECATED.md](./DEPRECATED.md).
* Compression should now be supported for all input and output files. Please [open an issue](https://github.com/nextstrain/augur/issues) if you find one that doesn't! [#1381][] (@victorlin)
* export v2: Add support to specify metadata columns to export without using them as colorings. This can be done with the `metadata_columns` property in the Auspice config JSON or via the `--metadata-columns` flag in the command line. [#1384][] (@joverlee521)

### Bug Fixes

* filter: In version 24.1.0, automatic conversion of boolean columns was accidentally removed. It has been restored with additional support for empty values evaluated as `None`. [#1410][] (@victorlin)
* filter: The order of rows in `--output-metadata` and `--output-strains` now reflects the order in the original `--metadata`. [#1294][] (@victorlin)
* filter, frequencies, refine: Performance improvements to reading the input metadata file. [#1294][] (@victorlin)
    * For filter, this comes with increased writing times for `--output-metadata` and `--output-strains`. However, net I/O speed still decreased during testing of this change.
* filter: Updated the help text of `--include` and `--include-where` to explicitly state that this can add strains that are missing an entry from `--sequences`. [#1389][] (@victorlin)
* filter: Fixed the summary messages to properly reflect force-inclusion of strains that are missing an entry from `--sequences`. [#1389][] (@victorlin)
* filter: Updated wording of summary messages. [#1389][] (@victorlin)
* Enforce UTF-8 encoding when reading and writing files. Improve error messages when a non-UTF-8 file is used. [#1381][] (@victorlin)

[#1294]: https://github.com/nextstrain/augur/pull/1294
[#1381]: https://github.com/nextstrain/augur/pull/1381
[#1384]: https://github.com/nextstrain/augur/pull/1384
[#1389]: https://github.com/nextstrain/augur/pull/1389
[#1410]: https://github.com/nextstrain/augur/pull/1410
[#1403]: https://github.com/nextstrain/augur/pull/1403

## 24.1.0 (30 January 2024)

### Features

* `augur.io.read_metadata`: A new optional `dtype` argument allows custom data types for all columns. Automatic type inference still happens by default, so this is not a breaking change. [#1252][] (@victorlin)
* `augur.io.read_vcf` has been removed and usage replaced with TreeTime's function of the same name which has improved validation of the VCF file. [#1366][] (@jameshadfield)

### Bug Fixes

* filter, frequencies, refine: Speed up reading of the metadata file. [#1252][] (@victorlin)
* traits: Previously, columns with only numeric values were treated as numerical data. These are now treated as categorical data for discrete trait analysis. [#1252][] (@victorlin)
* Support Biopython `≥1.82` by requiring bcbio-gff `≥0.7.1`. [#1400][] (@victorlin)

[#1252]: https://github.com/nextstrain/augur/pull/1252
[#1366]: https://github.com/nextstrain/augur/pull/1366
[#1400]: https://github.com/nextstrain/augur/pull/1400

## 24.0.0 (22 January 2024)

### Major Changes

* ancestral, translate: For VCF inputs please ensure you are using TreeTime 0.11.2 or later. A large number of bugfixes and improvements have been added in both Augur and TreeTime. [#1355][] and [TreeTime #263][] (@jameshadfield)
* ancestral, translate: GenBank files now require the (GFF mandatory) source feature to be present. [#1351][] (@jameshadfield)
* ancestral, translate: For GFF files, we extract the genome/sequence coordinates by inspecting the sequence-region pragma, region type and/or source type. This information is now required. [#1351][] (@jameshadfield)

### Features

* ancestral, translate: Improvements to VCF inputs / outputs. [#1355][] and [TreeTime #263][] (@jameshadfield)
    * Output VCF will better match the input VCF, including CHROM name and ploidy encoding.
    * VCF inputs now require `--vcf-reference-output`
    * AA sequences are now exported for the tree root
    * VCF writing is now 3 orders of magnitude faster (dataset dependent)
* ancestral, translate: A range of improvements to how we parse GFF and GenBank reference files. [#1351][] (@jameshadfield)
    * translate will now always export a 'nuc' annotation in the output JSON, allowing it to pass validation
    * Gene/CDS names of 'nuc' are now forbidden.
    * If a Gene/CDS in the GFF/GenBank file is unparsed we now print a warning.
* ancestral: For VCF alignments, a VCF output file is now only created when requested via `--output-vcf`. [#1344][] (@jameshadfield)
* ancestral: Improvements to command line arguments. [#1344][] (@jameshadfield)
     * Incompatible arguments are now checked, especially related to VCF vs FASTA inputs.
     * `--vcf-reference` and `--root-sequence` are now mutually exclusive.
* translate: Tree nodes are checked against the node-data JSON input to ensure sequences are present. [#1348][] (@jameshadfield)
* utils::load_features: This function may now raise `AugurError`. [#1351][] (@jameshadfield)
* export v2: Automatically minify large outputs. Use `--no-minify-json` to disable this default behavior. [#1352][] (@victorlin)
* Added a new file [DEPRECATED.md](./DEPRECATED.md) to document timelines and progress of deprecated features in the Augur CLI and Python API. [#1371][] (@victorlin)

### Bug Fixes

* ancestral, translate: Various fixes to VCF inputs / outputs. [#1355][] and [TreeTime #263][] (@jameshadfield)
    * Fix incorrect (but passing) tests
    * Fix case-sensitive sequence comparisons between the root and reference sequences.
    * Fix a bug where ambiguous alleles are not inferred (see [#1380][] for full details).
    * Fix a bug where positions with no sequence information were assigned a base because the mask was not being computed (see [#1382][] for full details).
    * More than one ALT allele is now correctly parsed
    * Mutations followed by an insertion are now parsed
    * Unchanged ref genotypes are now encoded as '0' rather than '.'
    * ALT alleles "*" are now valid (introduced in VCF spec 4.2, but observed in VCF 4.1 files)
    * Positions with no variation are no longer exported
* ancestral, translate: Fixes for JSON (non-VCF) inputs. [#1355][] (@jameshadfield)
    * The "reference" translations are now from the provided reference sequence, not from the root of the tree.  [#1355][] (@jameshadfield)
    * Fix a bug where positions with no sequence information were assigned a base because the mask was not applied (see [#1382][] for full details)
* ancestral, translate: Avoid incompatibilities with Biopython >=1.82. [#1374][], [#1387][] (@victorlin)
* ancestral, translate: Address Biopython deprecation warnings. [#1379][] (@victorlin)
* ancestral: Previously, the help text for `--genes` falsely claimed that it could accept a file. Now, it can truly claim that. [#1353][] (@victorlin)
* translate: The 'source' ID for GFF files is now ignored as a potential gene feature (it is still used for overall nuc coords). [#1348][] (@jameshadfield)
* translate: Improvements to command line arguments.  [#1348][] (@jameshadfield)
    * `--tree` and `--ancestral-sequences` are now required arguments.
    * separate VCF-only arguments into their own group
* translate: Fixes a bug in the parsing behaviour of GFF files whereby the presence of the `--genes` command line argument would change how we read individual GFF lines. Issue [#1349][], PR [#1351][] (@jameshadfield)
* If `TreeTimeError` is encountered Augur now exits with code 2 rather than 0. (This restores the original behaviour.) [#1367][] (@jameshadfield)
* Deprecate `read_strains` from `augur.utils` and add it to the public API under `augur.io`. [#1353][] (@victorlin)


[#1344]: https://github.com/nextstrain/augur/pull/1344
[#1348]: https://github.com/nextstrain/augur/pull/1348
[#1351]: https://github.com/nextstrain/augur/pull/1351
[#1349]: https://github.com/nextstrain/augur/issues/1349
[#1367]: https://github.com/nextstrain/augur/pull/1367
[#1371]: https://github.com/nextstrain/augur/pull/1371
[#1374]: https://github.com/nextstrain/augur/pull/1374
[#1379]: https://github.com/nextstrain/augur/pull/1379
[#1352]: https://github.com/nextstrain/augur/pull/1352
[#1353]: https://github.com/nextstrain/augur/pull/1353
[#1355]: https://github.com/nextstrain/augur/pull/1355
[#1380]: https://github.com/nextstrain/augur/issues/1380
[#1382]: https://github.com/nextstrain/augur/issues/1382
[#1387]: https://github.com/nextstrain/augur/pull/1387
[TreeTime #263]: https://github.com/neherlab/treetime/pull/263

## 23.1.1 (7 November 2023)

### Bug Fixes

* Fix Python 3.11 installation for Conda environments. [#1334][] (@victorlin)
* Bump `pyfastx` dependency to major versions 1 and 2. [#1335][] (@victorlin)

[#1334]: https://github.com/nextstrain/augur/issues/1334
[#1335]: https://github.com/nextstrain/augur/pull/1335

## 23.1.0 (22 September 2023)

### Features

* Support treetime 0.11.* [#1310][] (@corneliusroemer)
* export: Allow minimal export using only a (newick) tree in `augur export v2`. [#1299][] (@jameshadfield)
* A number of schema updates and improvements [#1299][] (@jameshadfield)
    * We now require all nodes to have `node_attrs` on them with one of `div` or `num_date` present
    * Some never-used properties are removed from the schemas, including a pattern for defining nucleotide INDELs which was never used by augur or auspice.
    * Tip label defaults are now settable within the auspice-config JSON
    * Empty colorings definitions are allowed (the tree will be grey in Auspice)

### Bug fixes

* ancestral: Export amino acid sequences inferred for the root node of the tree in the node data JSON output for compatibility with `augur translate` output. [#1317][] (@huddlej)

[#1299]: https://github.com/nextstrain/augur/pull/1299
[#1310]: https://github.com/nextstrain/augur/pull/1310
[#1317]: https://github.com/nextstrain/augur/pull/1317

## 23.0.0 (5 September 2023)

### Major Changes

* Drop support for Python 3.7. [#1296][] (@victorlin)

### Features

* export v2: Allow the root-sequence data to be included (inlined) in the main dataset JSON file, avoiding the need for a sidecar `_root-sequence.json` file. [#1295][] (@jameshadfield)

[#1295]: https://github.com/nextstrain/augur/pull/1295
[#1296]: https://github.com/nextstrain/augur/pull/1296

## 22.4.0 (29 August 2023)

### Features

* refine: Export covariance matrix and standard deviation for clock rate regression in the node data JSON output when these values are calculated by TreeTime. These new values appear in the `clock` data structure of the JSON output as `cov` and `rate_std` keys, respectively. [#1284][] (@huddlej)

### Bug fixes

* clades: Fix outputs for genes named `NA` (previously the value was replaced by `nan`). [#1293][] (@rneher)
* distance: Improve documentation by describing how gaps get treated as indels and how users can ignore specific characters in distance calculations. [#1285][] (@huddlej)
* Fix help output compatibility with non-Unicode streams. [#1290][] (@victorlin)

[#1284]: https://github.com/nextstrain/augur/pull/1284
[#1285]: https://github.com/nextstrain/augur/pull/1285
[#1290]: https://github.com/nextstrain/augur/pull/1290
[#1293]: https://github.com/nextstrain/augur/pull/1293

## 22.3.0 (14 August 2023)

### Features

* ancestral: add functionality to reconstruct ancestral amino acid sequences and add inferred mutations to the `node_data_json` with output equivalent to `augur translate`. `ancestral` now takes an annotation (`--annotation`), a list of genes (`--genes`), and a file name pattern for amino acid alignments (`--translations`). Mutations for each of these genes will be inferred and added to the output JSON to each node as a list at `['aa_muts'][gene]`. The annotations will be added to the `annotation` field in the output JSON. Inferred amino acids sequences can be saved with the new `--output-translations` argument. [#1258][] (@rneher, @huddlej)
* ancestral: add the ability to report mutations relative to a sequence other than the inferred root of the tree. This sequence can be specified via `--root-sequence` and difference between this sequence and the inferred root of the tree will be added as mutations to the root node for nucleotides and amino acids. All differences between the specified `root-sequence` and the inferred sequence of the root node of the tree will be added as mutations to the root node. This was previously already possible for `vcf` input via `--vcf-reference`. [#1258][] (@rneher)
* refine: add `mid_point` as rooting option to `refine`. [#1257][] (@rneher)

### Bug fixes

* filter: In version 22.2.0, `--query` would fail when the `.str` accessor was used on a column. This has been fixed. [#1277][] (@victorlin)

[#1257]: https://github.com/nextstrain/augur/pull/1257
[#1258]: https://github.com/nextstrain/augur/pull/1258
[#1277]: https://github.com/nextstrain/augur/issues/1277

## 22.2.0 (31 July 2023)

### Features

* Adds a new sub-command augur curate titlecase. The titlecase command is intended to apply titlecase to string fields in a metadata record (e.g. BRAINE-LE-COMTE, FRANCE -> Braine-le-Comte, France). Previously, this was available in the transform-string-fields script within the monkeypox repo.
 [#1197][] (@j23414 and @joverlee521)

[#1197]: https://github.com/nextstrain/augur/pull/1197

### Bug fixes

* export v2: Previously, when `strain` was not used as the metadata ID column, node attributes might have gone missing from the final Auspice JSON. This has been fixed. [#1260][], [#1262][] (@victorlin, @joverlee521)
* export v1: Added a deprecation warning for this command. [#1265][] (@victorlin)
* export v1: The recently introduced flag `--metadata-id-columns` did not work properly due to the same `export v2` bug that was fixed in this release. Instead of fixing it in `export v1`, drop the broken feature since this command is no longer being maintained. [#1265][] (@victorlin)
* filter: Expose internal Pandas errors from `--query` which may be useful to users. [#1267][] (@victorlin)
* filter: Previously, `--query` would fail when numerical comparisons were used on columns with missing values. This has been fixed. [#1269][] (@victorlin)

[#1260]: https://github.com/nextstrain/augur/issues/1260
[#1262]: https://github.com/nextstrain/augur/issues/1262
[#1265]: https://github.com/nextstrain/augur/pull/1265
[#1267]: https://github.com/nextstrain/augur/pull/1267
[#1269]: https://github.com/nextstrain/augur/issues/1269

## 22.1.0 (10 July 2023)

### Features

* export, frequencies, refine, traits: Add a new flag `--metadata-id-columns` to customize the possible metadata ID columns. Previously, this was only available in `augur filter`. [#1240][] (@victorlin)
* Add new sub-subcommand augur curate format-dates. The format-dates command is intended to be used to format date fields to ISO 8601 date format (YYYY-MM-DD), where incomplete dates are masked with `XX` (e.g. 2023 -> 2023-XX-XX). [#1146][] (@joverlee521)

### Bug fixes

* parse: Fix a bug where `--fix-dates` was always applied, with a default of `--fix-dates=monthfirst`. Now, running without `--fix-dates` will leave dates as-is. [#1247][] (@victorlin)
* `augur.io.open_file`: Previously, the docs described a type restriction on `path_or_buffer` but it was not enforced. It has been updated to allow all I/O classes, and is enforced at run-time. [#1250][] (@victorlin)
* filter: Fix a bug where data files consisting of only numerical strain names would not work when both `--metadata` and `--sequences` are passed. [#1256][] (@victorlin)

[#1146]: https://github.com/nextstrain/augur/pull/1146
[#1240]: https://github.com/nextstrain/augur/pull/1240
[#1247]: https://github.com/nextstrain/augur/issues/1247
[#1250]: https://github.com/nextstrain/augur/pull/1250
[#1256]: https://github.com/nextstrain/augur/pull/1256

## 22.0.3 (14 June 2023)

### Bug fixes

* utils: Serialize pandas Series in `write_json`. [#1213][] (@victorlin)

[#1213]: https://github.com/nextstrain/augur/pull/1213

## 22.0.2 (26 May 2023)

### Bug fixes

* CI: Add a Github action to test augur on 8 Nextstrain pathogen workflows using example data. [#1217][] (@corneliusroemer)
* parse: Denote required arguments including `--fields`, `--output-sequences`, and `--output-metadata`. [#1228][] (@huddlej)
* Fix export of the `strand` attribute of gene annotations. Previously, features on the negative strand were not annotated as such since the code assumed that the `strand` attribute was boolean instead of `[-1, +1]`. [#1211] @rneher and @j23414.
* augur.io.read_metadata: explicitly set `date` column as `string` type to prevent year only dates from being inferred as integers. [#1235][] (@joverlee521)

[#1211]: https://github.com/nextstrain/augur/pull/1211
[#1217]: https://github.com/nextstrain/augur/pull/1217
[#1228]: https://github.com/nextstrain/augur/pull/1228
[#1235]: https://github.com/nextstrain/augur/pull/1235

## 22.0.1 (16 May 2023)

### Bug fixes

* export: No longer export duplicate entries in the colorings array, a bug which has been present in Augur since at least v12 [#719][]. [#1218][] (@jameshadfield)
* export: In version 22.0.0, some configurations of export may have resulted in the clade coloring appearing last in the Auspice dropdown rather than first. This is now fixed. [#1218] (@jameshadfield)
* export: In version 22.0.0, validation of `augur.utils.read_node_data` was changed to error when a node data JSON did not contain any actual data. This causes export to error when an empty node data JSON is passed, as for example in ncov's pathogen-ci. This is now fixed by warning instead. The bug was originally introduced in PR [#728][]. [#1214][] (@corneliusroemer)

[#719]: https://github.com/nextstrain/augur/issues/719
[#1214]: https://github.com/nextstrain/augur/pull/1214
[#1218]: https://github.com/nextstrain/augur/pull/1218

## 22.0.0 (9 May 2023)

### Major Changes

* export, filter, frequencies, refine, traits: From versions 10.0.0 through 21.1.0, arbitrary delimiters for `--metadata` were supported due to internal implementation differences from the advertised CSV and TSV support. Starting with this version, non-CSV/TSV files will no longer be supported by default. To adjust for this breaking change, specify custom delimiters with the new `--metadata-delimiters` flag. [#1196][] (@victorlin)
* `augur.io.read_metadata`: Previously, this supported any arbitrary delimiters for the metadata. Now, it only supports a list of possible delimiters represented by the new `delimiters` keyword argument, which defaults to `,` and `\t`. [#812][] (@victorlin)
* refine: The seeding method for `--seed` has been updated. This affects usages that rely on the reproducibility of outputs with the same `--seed` value prior to this version. Outputs from this version onwards should be reproducible until the next implementation change, which we don't expect to happen any time soon. [#1207][] (@rneher)

### Features

* Constrain `bcbio-gff` to >=0.7.0 and allow `Biopython` >=1.81 again. We had to introduce the `Biopython` constraint in v21.0.1 (see [#1152][]) due to `bcbio-gff` <0.7.0 relying on the removed `Biopython` feature `UnknownSeq`. [#1178][] (@corneliusroemer)
* `augur.io.read_metadata` (used by export, filter, frequencies, refine, and traits): Previously, this used the Python parser engine for [`pandas.read_csv()`][]. Updated to use the C engine for faster reading of metadata. [#812][] (@victorlin)
* curate: Allow custom metadata delimiters with the new `--metadata-delimiters` flag. [#1196][] (@victorlin)
* Bump the default recursion limit to 10,000. Users can continue to override this limit with the environment variable `AUGUR_RECURSION_LIMIT`. [#1200][] (@joverlee521)
* clades, export v2: Clade labels + coloring keys are now definable via arguments to augur clades allowing pipelines to use multiple invocations of augur clades resulting in multiple sets of colors and branch labels. How labels are stored in the (intermediate) node-data JSON files has changed. This should be fully backwards compatible for pipelines using augur commands, however custom scripts may need updating. PR [#728][] (@jameshadfield)
* refine: add flag `--max-iter` to control the maximal number of iterations TreeTime uses to infer time trees. This was previously hard-coded to 2, which is now the default. [#1203][] (@rneher)
* refine: add flags `--greedy-resolve` and `--stochastic-resolve` to customize polytomy resolution. [#1203][], [#1207][] (@rneher)
  * `--greedy-resolve`: resolve polytomies by greedily minimizing tree length (default behavior, unchanged).
  * `--stochastic-resolve`: resolve polytomies as random coalescent trees.
  * These are mutually exclusive with the pre-existing `--keep-polytomies` flag.

### Bug fixes

* filter, frequencies, refine, parse: Previously, ambiguous dates in the future had a limit of today's date imposed on the upper value but not the lower value. It is now imposed on the lower value as well. [#1171][] (@victorlin)
* refine: `--year-bounds` was ignored in versions 9.0.0 through 20.0.0. It now works. [#1136][] (@victorlin)
* tree: Input alignment filenames which do not end in `.fasta` are now properly handled when using IQ-TREE.  Previously their contents were overwritten first by `augur tree` itself (resulting in truncation) and then by the log output of IQ-TREE (resulting in an error).  Thanks to Jon Bråte for reporting this bug. [#1206][] (@tsibley)
* clades: A number of small bug fixes, improvements to documentation, tests and improved error detection. [#1199][] (@jameshadfield)

[#728]: https://github.com/nextstrain/augur/pull/728
[#812]: https://github.com/nextstrain/augur/pull/812
[#1136]: https://github.com/nextstrain/augur/issues/1136
[#1152]: https://github.com/nextstrain/augur/pull/1152
[#1171]: https://github.com/nextstrain/augur/issues/1171
[#1178]: https://github.com/nextstrain/augur/pull/1178
[#1196]: https://github.com/nextstrain/augur/pull/1196
[#1199]: https://github.com/nextstrain/augur/pull/1199
[#1200]: https://github.com/nextstrain/augur/pull/1200
[#1203]: https://github.com/nextstrain/augur/pull/1203
[#1206]: https://github.com/nextstrain/augur/pull/1206
[#1207]: https://github.com/nextstrain/augur/pull/1207
[`pandas.read_csv()`]: https://pandas.pydata.org/pandas-docs/version/1.5/reference/api/pandas.read_csv.html

## 21.1.0 (14 March 2023)

### Features

* filter: Add `--empty-output-reporting={error,warn,silent}` option to allow filter to produce empty outputs without raising an error. The default behavior is still to raise an error when filter produces an empty output, so users will have to explicitly pass the "warn" or "silent" value to bypass the error. [#1175][] (@joverlee521)

### Bug fixes

* translate: Fix error handling when features cannot be read from reference sequence file. [#1168][] (@victorlin)
* translate: Remove an unnecessary check which allowed for inaccurate error messages to be shown. [#1169][] (@victorlin)
* frequencies: Previously, monthly pivot points calculated from the end of a month may have been shifted by 1-3 days. This is now fixed. [#1150][] (@victorlin)
* docs: Fix minor formatting issues. [#1095][] (@victorlin)
* Update development status on PyPI from "3 - Alpha" to "5 - Production/Stable". This should have been done since the beginning of this changelog, but now it is official. [#1160][] (@corneliusroemer)

[#1095]: https://github.com/nextstrain/augur/pull/1095
[#1150]: https://github.com/nextstrain/augur/pull/1150
[#1160]: https://github.com/nextstrain/augur/pull/1160
[#1168]: https://github.com/nextstrain/augur/pull/1168
[#1169]: https://github.com/nextstrain/augur/pull/1169
[#1175]: https://github.com/nextstrain/augur/pull/1175

## 21.0.1 (17 February 2023)

### Bug fixes

* Constrain Biopython version to <=1.80 so that `augur translate` is not broken by a deprecation of `UnknownSeq` in 1.81. When running `augur translate` with Biopython 1.81, the user will receive an error starting with `ERROR: Package BCBio.GFF not found!` and ending with `TypeError: object of type 'NoneType' has no len()`. [#1152][] (@corneliusroemer)

[#1152]: https://github.com/nextstrain/augur/pull/1152

## 21.0.0 (7 February 2023)

### Major Changes

* measurements export: Supports exporting multiple thresholds per collection via the measurements config and the `--thresholds` option. This change is backwards compatible with previous uses of the `--threshold` option. However, due to the updates to the JSON schema, users will need to update to Auspice v2.43.0 for thresholds to be displayed properly in the measurements panel. [#1148][] (@joverlee521)

### Features

* export v2: Add `--validation-mode={error,warn,skip}` option for more nuanced control of validation.  The new "warn" mode performs validation and emits messages about potential problems, but it does not cause the export command to fail even if there are problems. [#1135][] (@tsibley)

### Bug Fixes

* filter, frequencies, refine, parse: Properly handle invalid date errors and output the bad date. [#1140][] (@victorlin)
* export, validate: Validation errors are now much more human-readable and actually pinpoint the problems. [#1134][] (@tsibley)

[#1134]: https://github.com/nextstrain/augur/pull/1134
[#1135]: https://github.com/nextstrain/augur/pull/1135
[#1140]: https://github.com/nextstrain/augur/pull/1140
[#1148]: https://github.com/nextstrain/augur/pull/1148

## 20.0.0 (20 January 2023)

### Major Changes

* frequencies: Changes the logic for calculating the time points when frequencies are estimated to ensure that the user-provided "end date" is always included. This change in the behavior of the frequencies command fixes a bug where large intervals between time points (e.g., 3 months) could cause recent data to be omitted from frequency calculations. See the pull request for more details included the scientific implications of this bug. [#1121][] (@huddlej)

[#1121]: https://github.com/nextstrain/augur/pull/1121

## 19.3.0 (19 January 2023)

### Features

* titers: Support parsing of thresholded values (e.g., "<80" or ">2560"). [#1118][] (@huddlej)
* tree: Support bootstrapped trees generated with RAxML via user-provided `--tree-builder-args`. [#1127][] (@tsibley)

### Bug Fixes

* utils: Serialize common numpy data types in `write_json`. [#1119][] (@victorlin)
* filter: Standardize exit codes from internal error handling. [#931][] (@victorlin)
* tree: Suppress the `Cannot specify --substitution-model unless using IQTree` warning when `--substitution-model` is left at its default. [#1127][] (@tsibley)
* tree: Print the underlying error message when tree building fails. [#1127][] (@tsibley)
* Previously, `numpy` and `scipy` were installed as dependencies of dependencies. Mark them as direct dependencies since they are used directly within Augur. [#1120][] (@victorlin)

[#931]: https://github.com/nextstrain/augur/pull/931
[#1118]: https://github.com/nextstrain/augur/pull/1118
[#1119]: https://github.com/nextstrain/augur/pull/1119
[#1120]: https://github.com/nextstrain/augur/pull/1120
[#1127]: https://github.com/nextstrain/augur/pull/1127

## 19.2.0 (19 December 2022)

### Features

* titers: Allow users to specify a custom prefix for attributes in the JSON output (e.g., `cTiter` can be changed to `custom_prefix_cTiter`). [#1106][] (@huddlej)

[#1106]: https://github.com/nextstrain/augur/pull/1106

## 19.1.0 (14 December 2022)

### Features

* io: Add `open_file` and `write_sequences` to the Python Pubic API. [#1114][] (@joverlee521)

[#1114]: https://github.com/nextstrain/augur/pull/1114

## 19.0.0 (13 December 2022)

### Major Changes

* io: Only `read_metadata` and `read_sequences` are available as part of the Python Public API. Other Python API functions of the `augur.io` module are no longer directly available. This is a breaking change, although we suspect few users to be impacted. If you still need to use other imports in your scripts, they can be imported from the [Developer API](https://docs.nextstrain.org/projects/augur/en/stable/api/developer/index.html) but note that they are no longer part of the [Public API](https://docs.nextstrain.org/projects/augur/en/stable/api/public/index.html). [#1087][] (@victorlin)

### Bug Fixes

* docs: Update the API documentation to reflect the latest state of things in the codebase. [#1087][] (@victorlin)
* Fix support for Biopython version 1.80 which deprecated `Bio.Seq.Seq.ungap()`. [#1102][] (@victorlin)
* export v2: Fixed a bug where colorings for zero values via `--colors` would not get applied to the exported Auspice JSON. [#1100][] (@joverlee521)
* curate: Fixed a bug where metadata TSVs failed to parse if data within a column included comma separated values [#1110][] (@joverlee521)

[#1087]: https://github.com/nextstrain/augur/pull/1087
[#1100]: https://github.com/nextstrain/augur/pull/1100
[#1102]: https://github.com/nextstrain/augur/pull/1102
[#1110]: https://github.com/nextstrain/augur/pull/1110

## 18.2.0 (15 November 2022)

### Features

* Add the curate subcommand with two sub-subcommands, passthru and normalize-strings. The curate subcommand is intended to be a suite of commands to help users with data curation prior to running Nextstrain analyses. We will continue to add more subcommands as we identify other common data curation tasks. Please see the [usage docs](https://docs.nextstrain.org/projects/augur/en/stable/usage/cli/curate) for details. [#1039][] (@joverlee521)

[#1039]: https://github.com/nextstrain/augur/pull/1039

## 18.1.2 (1 November 2022)

### Bug Fixes

* traits: Fix trait inference when tips have missing values. [#1081][] (@huddlej)

[#1081]: https://github.com/nextstrain/augur/pull/1081

## 18.1.1 (1 November 2022)

### Bug Fixes

* filter: Fixed a bug where `--group-by week` would fail when all samples in a chunk have been dropped due to ambiguous dates. [#1080][] (@victorlin)

[#1080]: https://github.com/nextstrain/augur/pull/1080

## 18.1.0 (26 October 2022)

### Features

* filter: Add support to group by ISO week (`--group-by week`) during subsampling. [#1067][] (@victorlin)

### Bug Fixes

* filter: Fixed unintended behavior in which grouping by `day` would "work" when used with `month` and/or `year`. Updated so it will be ignored. [#1070][] (@victorlin)
* filter: Fixed unintended behavior in which grouping by `month` with ambiguous years would "work". Updated so date ambiguity is checked properly for all generated columns. [#1072][] (@victorlin)

[#1067]: https://github.com/nextstrain/augur/pull/1067
[#1070]: https://github.com/nextstrain/augur/pull/1070
[#1072]: https://github.com/nextstrain/augur/pull/1072

## 18.0.0 (21 September 2022)

### Major Changes

* export: The `--node-data` option may now be given multiple times to provide additional `.json` files.  Previously, subsequent occurrences of the option overrode prior occurrences.  This is a **breaking change**, although we expect few usages to be impacted.  Each occurrence of the option may still specify multiple files at a time. [#1010][] (@tsibley)

### Bug Fixes


* refine: 17.1.0 updated TreeTime to version 0.9.2 and introduced the `refine` flag `--use-fft`. This makes previously costly marginal date inference cheaper. This update adjusts when `refine` runs marginal date inference during its iterative optimization. Without the `use-fft` flag, it will now behave as it did before 17.1.0 (marginal inference only during final iterations). With the `--use-fft` flag, marginal date inference will be used at every step during the iteration if refine is run with `--date-inference marginal` [#1034][]. (@rneher)
* tree: When using IQtree as tre builder, `--nthreads` now sets the maximum number of threads (IQtree argument `-ntmax`). The actual number of threads to use can be specified by the user through the tree-builder-arg `-nt` which defaults to `-nt AUTO`, causing IQtree to automatically chose the best number of threads to use [#1042][] (@corneliusroemer)
* Make cvxopt as a required dependency, since it is required for titer models to work [#1035][]. (@victorlin)
* filter: Fix compatibility with Pandas 1.5.0 which could cause an unexpected `AttributeError` with an invalid `--query` given to `augur filter`. [#1050][] (@tsibley)
* refine: Add `--verbosity` argument that is passed down to TreeTime to facilitate monitoring and debugging. [#1033][] (@anna-parker)
* Improve handling of errors from TreeTime. [#1033][] (@anna-parker)

[#1010]: https://github.com/nextstrain/augur/pull/1010
[#1033]: https://github.com/nextstrain/augur/pull/1033
[#1034]: https://github.com/nextstrain/augur/pull/1034
[#1035]: https://github.com/nextstrain/augur/pull/1035
[#1042]: https://github.com/nextstrain/augur/pull/1042
[#1050]: https://github.com/nextstrain/augur/pull/1050


## 17.1.0 (19 August 2022)

### Features

* refine: Upgrade TreeTime from 0.8.6 to >= 0.9.2 which enables a speedup of timetree inference in marginal mode due to the use of Fast Fourier Transforms [#1018][]. (@rneher and @anna-parker). Use the `refine` flag `--use-fft` to use this feature.

### Bug Fixes

* refine, export v1: Use pandas.DataFrame.at instead of .loc for single values [#979][]. (@victorlin)
* refine: Gracefully handle all exceptions from TreeTime [#1023][]. (@anna-parker)
* refine: Document branch length units `treetime` expects [#1024][]. (@anna-parker)
* dates: Raise an error when metadata to `get_numerical_dates()` is not a pandas DataFrame [#1026][]. (@victorlin)

[#979]: https://github.com/nextstrain/augur/pull/979
[#1018]: https://github.com/nextstrain/augur/pull/1018
[#1023]: https://github.com/nextstrain/augur/pull/1023
[#1024]: https://github.com/nextstrain/augur/pull/1024
[#1026]: https://github.com/nextstrain/augur/pull/1026

## 17.0.0 (9 August 2022)

### Major Changes

* Moved the following modules to subpackages [#1002][]. (@joverlee521)
  These are technically breaking changes for the API, but they do not change the Augur CLI commands.
    * `import.py` -> `import_/__init__.py`
    * `import_beast.py` -> `import_/beast.py`
    * `measurements.py` -> `measurements/__init__.py` + `measurements/concat.py` + `measurements/export.py`
* Move the following internal functions/classes [#1002][]. (@joverlee521)
    * `augur.add_default_command` -> `argparse_.add_default_command`
    * `utils.HideAsFalseAction` -> `argparse_.HideAsFalseAction`
* Subcommands must include a `register_parser` function to add their own parser instead of a `register_arguments` function [#1002][]. (@joverlee521)
* utils: Remove internal function `utils.read_metadata()` [#978][]. (@victorlin)
    * Use `io.read_metadata()` going forwards.
    * To switch to using metadata as a pandas DataFrame (recommended):
        * Iterate through strains: `metadata.items()` -> `metadata.iterrows()`
        * Check strain presence: `strain in metadata` -> `strain in metadata.index`
        * Check field presence: `field in metadata[strain]` -> `field in metadata.columns`
        * Get metadata for a strain: `metadata[strain]` -> `metadata.loc[strain]`
        * Get field for a strain: `metadata[strain][field]` -> `metadata.at[strain, field]`
    * To keep using metadata in a dictionary:
        ```py
        metadata = read_metadata(args.metadata)
        metadata.insert(0, "strain", metadata.index.values)
        columns = metadata.columns
        metadata = metadata.to_dict(orient="index")
        ```

### Features

* export: `--skip-validation` now also skips version compatibility checks [#902][]. (@corneliusroemer)
* filter: Report names of duplicate strains found during metadata parsing [#1008][] (@huddlej)
* translate: Add support for Nextclade gene map GFFs [#1017][] (@huddlej)

### Bug Fixes

* filter: Rename internal force inclusion filtering functions [#1006][] (@victorlin)

[#902]: https://github.com/nextstrain/augur/pull/902
[#978]: https://github.com/nextstrain/augur/pull/978
[#1002]: https://github.com/nextstrain/augur/pull/1002
[#1006]: https://github.com/nextstrain/augur/pull/1006
[#1008]: https://github.com/nextstrain/augur/pull/1008
[#1017]: https://github.com/nextstrain/augur/pull/1017

## 16.0.3 (6 July 2022)

### Bug Fixes

* filter: Move `register_arguments` to the top of the module for better readability [#995][]. (@victorlin)
* filter: Fix a regression [introduced in 16.0.2](https://github.com/nextstrain/augur/commit/4859b5d70e77cc9a0bb99e741fefb29952058b71) that caused grouping with subsampled max sequences and force-included strains to fail in a data-specific way [#1000][]. (@huddlej)

[#995]: https://github.com/nextstrain/augur/pull/995
[#1000]: https://github.com/nextstrain/augur/pull/1000

## 16.0.2 (30 June 2022)

### Bug Fixes

* The entropy panel was unavailable if mutations were not translated [#881][]. This has been fixed by creating an additional `annotations` block in `augur ancestral` containing (nucleotide) genome annotations in the node-data [#961][] (@jameshadfield)
* ancestral: WARNINGs to stdout have been updated to print to stderr [#961][] (@jameshadfield)
* filter: Explicitly drop date/year/month columns from metadata during grouping. [#967][] (@victorlin)
    * This fixes a bug [#871][] where `augur filter` would crash with a cryptic `ValueError` if `year` and/or `month` is a custom column in the input metadata and also included in `--group-by`.
* filter: Fix duplicates that may appear in metadata when using `--include`/`--include-where` with subsampling [#986][] (@victorlin)

[#881]: https://github.com/nextstrain/augur/issues/881
[#961]: https://github.com/nextstrain/augur/pull/961
[#967]: https://github.com/nextstrain/augur/pull/967
[#871]: https://github.com/nextstrain/augur/issues/871
[#986]: https://github.com/nextstrain/augur/pull/986

## 16.0.1 (21 June 2022)

### Bug Fixes

* filter: Handle errors from `filter_by_query` [#942][] (@victorlin)
* translate: output nuc annotation when reading from gff3 gene map [#976][] (@corneliusroemer)
* CI: Remove step for selecting PyPI instance [#974][] (@victorlin)
* CI: Add token to use GitHub CLI [#958][] (@victorlin)

[#942]: https://github.com/nextstrain/augur/pull/942
[#976]: https://github.com/nextstrain/augur/pull/976
[#974]: https://github.com/nextstrain/augur/pull/974
[#958]: https://github.com/nextstrain/augur/pull/958

## 16.0.0 (16 June 2022)

### Major Changes

* filter: Error when any group-by column is not found [#933][] (@victorlin)
    * Check your workflows for any new errors that may arise from this.
* parse: Error on duplicates instead of silently passing [#918][] (@victorlin)
    * Check your workflows for any new errors that may arise from this.
* utils: Remove `utils.myopen()` [#926][] (@victorlin)
    * Use `io.open_file()` going forwards.
* Moved the following internal functions [#929][], [#923][] (@victorlin):
    * `utils.read_vcf` -> `io.read_vcf`
    * `utils.run_shell_command` -> `io.run_shell_command`
    * `utils.shquote` -> `io.shquote`
    * `utils.ambiguous_date_to_date_range` -> `dates.ambiguous_date_to_date_range`
    * `utils.is_date_ambiguous` -> `dates.is_date_ambiguous`
    * `utils.get_numerical_date_from_value` -> `dates.get_numerical_date_from_value`
    * `utils.get_numerical_dates` -> `dates.get_numerical_dates`
        * Drop support for dict type as the first parameter [#934][]
    * `filter.write_vcf` -> `io.write_vcf`

### Features

* Add the measurements subcommand with two sub-subcommands, export and concat [#879][] (@joverlee521)
* filter: Report min and max date separately [#930][] (@victorlin)
* export v2: Allow the color scale type to be temporal [#969][] (@jameshadfield)
* Handle `FileNotFoundError` and unexpected exceptions gracefully [#914][] (@victorlin)

### Bug Fixes

* filter: Properly handle error on duplicates [#918][] (@victorlin)
* filter: Reorganize Cram test files [#943][] (@victorlin)
* filter: Reword comment on vcftools [#924][] (@victorlin)
* io: Split io.py into smaller files under new io/ [#949][] (@victorlin)
* io: Add tests for `io.open_file()` [#926][] (@victorlin)
* Move AugurError to new errors.py, replace RuntimeError [#921][] (@victorlin)
* Remove internal usage of `utils.read_metadata()` [#934][], [#972][] (@victorlin)
* schemas: Add missing display_default properties for Auspice config v2 [#916][] (@tsibley)
* CI: Split codecov into separate job, combine coverage of all matrix jobs [#968][] (@tsibley)
* CI: Temporarily disable failing test [#962][] (@victorlin)
* CI: pip install without editable mode [#956][] (@victorlin)
* CI: Include functional tests in code coverage [#899][] (@huddlej)
* CI: Move --quiet flag to accommodate snakemake=7.7.0 behavior [#927][] (@victorlin)
* CI: Move docker rebuild step to release workflow [#912][] (@victorlin)
* Update release process [#913][] (@victorlin)

[#968]: https://github.com/nextstrain/augur/pull/968
[#949]: https://github.com/nextstrain/augur/pull/949
[#879]: https://github.com/nextstrain/augur/pull/879
[#962]: https://github.com/nextstrain/augur/pull/962
[#956]: https://github.com/nextstrain/augur/pull/956
[#921]: https://github.com/nextstrain/augur/pull/921
[#943]: https://github.com/nextstrain/augur/pull/943
[#929]: https://github.com/nextstrain/augur/pull/929
[#934]: https://github.com/nextstrain/augur/pull/934
[#930]: https://github.com/nextstrain/augur/pull/930
[#926]: https://github.com/nextstrain/augur/pull/926
[#933]: https://github.com/nextstrain/augur/pull/933
[#899]: https://github.com/nextstrain/augur/pull/899
[#923]: https://github.com/nextstrain/augur/pull/923
[#914]: https://github.com/nextstrain/augur/pull/914
[#913]: https://github.com/nextstrain/augur/pull/913
[#927]: https://github.com/nextstrain/augur/pull/927
[#924]: https://github.com/nextstrain/augur/pull/924
[#916]: https://github.com/nextstrain/augur/pull/916
[#918]: https://github.com/nextstrain/augur/pull/918
[#912]: https://github.com/nextstrain/augur/pull/912
[#969]: https://github.com/nextstrain/augur/pull/969
[#972]: https://github.com/nextstrain/augur/pull/972

## 15.0.2 (5 May 2022)

### Bug Fixes

* docs: Fix API documentation rendering and add page for `io` module [#896][] (@joverlee521)
* CI: Use GitHub Actions for release process [#904][] (@victorlin)
* utils: Fix branch length annotations in `json_to_tree` function [#908][] (@huddlej)
* export v2: Use io.read_metadata during export, fixing a bug caused when the user's input metadata does not have any valid strain id columns [#909][] (@huddlej)
* CI: Call new GitHub Actions workflow to rebuild images [#910][] (@victorlin)

[#910]: https://github.com/nextstrain/augur/pull/910
[#909]: https://github.com/nextstrain/augur/pull/909
[#908]: https://github.com/nextstrain/augur/pull/908
[#904]: https://github.com/nextstrain/augur/pull/904
[#896]: https://github.com/nextstrain/augur/pull/896

## 15.0.1 (25 April 2022)

### Bug Fixes

* traits: Fix crash when inferring traits with a single value. [#891][] (@huddlej)
* index: Properly handle missing inputs. [#900][] (@huddlej)
* index: Use standard UNIX-style line endings for output. [#900][] (@huddlej)

[#891]: https://github.com/nextstrain/augur/pull/891
[#900]: https://github.com/nextstrain/augur/pull/900

## 15.0.0 (15 April 2022)

### Major Changes

* export: Move extensions block to meta [#888][] (@corneliusroemer)

    Note: this is technically a breaking change, but the misplaced extensions block was added in version 14.1.0 and intended for internal use by Nextclade. We don't expect any users to be impacted by this.

### Features

* filter: Support relative dates for `--min-date` and `--max-date` [#740][] (@benjaminotter, @victorlin)
* frequencies: Support relative dates for `--min-date` and `--max-date` [#890][] (@victorlin, @huddlej)
* docs: Add page describing how to use your own tree builder [#884][] (@jameshadfield, @joverlee521)

### Bug Fixes

* parse: Fix `--fix-dates` [#878][] (@huddlej)
* tree: Fix internal comment on excluded sites [#880][] (@jameshadfield)

[#888]: https://github.com/nextstrain/augur/pull/888
[#740]: https://github.com/nextstrain/augur/pull/740
[#890]: https://github.com/nextstrain/augur/pull/890
[#884]: https://github.com/nextstrain/augur/pull/884
[#878]: https://github.com/nextstrain/augur/pull/878
[#880]: https://github.com/nextstrain/augur/pull/880

## 14.1.0 (31 March 2022)

### Features

* schemas: Extend export v2 schema to support an array of trees [#851][] (@tsibley)
* schemas: Add JSON schemas for our root-sequence and tip-frequencies sidecars [#852][] (@tsibley)
* schemas: Add JSON schema for measurements sidecar [#859][] (@joverlee521)
* filter: Send warnings to stderr to be consistent with other warnings [#862][] (@victorlin)
* export: Allow an extensions block in auspice config & dataset JSONs [#865][] (@jameshadfield)
* export: Allow skipping of input/output schema validation [#865][] (@jameshadfield)
* export:  Order keys in dataset for easier reading [#868][] (@jameshadfield)

### Bug Fixes

* parse: Fix typo in internal variable name [#850][] (@emmahodcroft)

[#850]: https://github.com/nextstrain/augur/pull/850
[#851]: https://github.com/nextstrain/augur/pull/851
[#852]: https://github.com/nextstrain/augur/pull/852
[#859]: https://github.com/nextstrain/augur/pull/859
[#862]: https://github.com/nextstrain/augur/pull/862
[#865]: https://github.com/nextstrain/augur/pull/865
[#868]: https://github.com/nextstrain/augur/pull/868

## 14.0.0 (8 February 2022)

### Major Changes

* Drop support for Python 3.6, add support for 3.9 and 3.10 [#822][] (@victorlin)

### Features

* refine: Enable bootstrap support by passing confidence values through to Auspice JSONs from `augur refine` node data JSONs [#839][] (@huddlej)
* tree: Allow users to override default tree builder arguments with a new `--override-default-args` flag [#839][] (@huddlej)
* clades: Allow descendant clades to be defined by explicitly inheriting from ancestral clade names [#846][] (@corneliusroemer)

### Bug Fixes

* tree: Fix segmentation fault that can occur when user-provided tree builder args conflict hardcoded defaults for IQ-TREE's. The new `--override-default-args` flag allows users to override the defaults that conflict with their values. [#839][] (@huddlej)
* filter/utils: fix year-only and numeric date handling [#841][] (@victorlin)
* CI: test earliest supported Biopython versions in matrix, remove redundant installs [#843][] (@victorlin)

[#822]: https://github.com/nextstrain/augur/pull/822
[#841]: https://github.com/nextstrain/augur/pull/841
[#843]: https://github.com/nextstrain/augur/pull/843
[#846]: https://github.com/nextstrain/augur/pull/846
[#839]: https://github.com/nextstrain/augur/pull/839

## 13.1.2 (28 January 2022)

### Features

### Bug Fixes

* filter: Set default filter priority to negative infinity instead of 0. [#835][] (@huddlej)
* filter: Cast strains to list when subsetting metadata to prevent pandas FutureWarning. [#833][] (@victorlin)
* export v2: Allow periods (.) in accession numbers. [#832][] (@tsibley)

[#835]: https://github.com/nextstrain/augur/pull/835
[#833]: https://github.com/nextstrain/augur/pull/833
[#832]: https://github.com/nextstrain/augur/pull/832

## 13.1.1 (21 January 2022)

### Bug Fixes

* dependencies: Relax upper bounds on development dependencies. [#821][] and [#830][] (@victorlin, @huddlej)
* refine: Inform users when their requested root sequence is not present in the alignment. [#816][] (@corneliusroemer)
* utils: Fix crash with `read_metadata` when numexpr is installed. [#826][] (@victorlin)

[#816]: https://github.com/nextstrain/augur/pull/816
[#821]: https://github.com/nextstrain/augur/pull/821
[#826]: https://github.com/nextstrain/augur/pull/826
[#830]: https://github.com/nextstrain/augur/pull/830

## 13.1.0 (10 December 2021)

### Features

* schemas: Add "$id" key to Auspice config schemas so we have a way of referring to these. [#806][] (@tsibley)

### Bug Fixes

* filter: Fix groupby with incomplete dates. [#808][] (@victorlin)

[#806]: https://github.com/nextstrain/augur/pull/806
[#808]: https://github.com/nextstrain/augur/pull/808

## 13.0.4 (8 December 2021)

### Bug Fixes

* dependencies: Replace deprecated mutable sequence interface for BioPython. [#788][] (@Carlosbogo)
* dependencies: Fix backward compatibility with BioPython. [#801][] (@huddlej)
* data: Add latitude and longitude details for "Reunion". [#791][] (@corneliusroemer)
* filter: Use pandas functions to determine subsample groups. [#794][] and [#797][] (@victorlin)
* filter: Add clarity to help message and output of probabilistic sampling. [#792][] (@victorlin)

[#788]: https://github.com/nextstrain/augur/pull/788
[#791]: https://github.com/nextstrain/augur/pull/791
[#792]: https://github.com/nextstrain/augur/pull/792
[#794]: https://github.com/nextstrain/augur/pull/794
[#797]: https://github.com/nextstrain/augur/pull/797
[#801]: https://github.com/nextstrain/augur/pull/801

## 13.0.3 (19 November 2021)

### Bug Fixes

* tree: Handle compressed alignment when excluding sites. [#786][] (@huddlej)
* docs: Fix typos ([ce0834c][]) and clarify exclude sites inputs ([5ad1574][]). (@corneliusroemer)

[#786]: https://github.com/nextstrain/augur/pull/786
[ce0834c]: https://github.com/nextstrain/augur/commit/ce0834c476abc9ee99785fa930608218b7d78990
[5ad1574]: https://github.com/nextstrain/augur/commit/5ad157485015623883c6b637d247459f906b63cb

## 13.0.2 (12 October 2021)

### Bug Fixes

* dependencies: Support latest versions of BioPython. [#777][] (@huddlej)
* tree: Allow users to specify arbitrary IQ-TREE models. [#776][] (@huddlej)

[#776]: https://github.com/nextstrain/augur/pull/776
[#777]: https://github.com/nextstrain/augur/pull/777

## 13.0.1 (1 October 2021)

### Bug Fixes

* docs: Fix broken link to latitude/longitude documentation. [#766][] (@victorlin)
* filter: Fix reproducibility of subsampling by using the user-defined random seed in all random function calls and by ordering strain sets as lists prior to adding strains to group-by priority queues. [#772][] (@huddlej)

[#766]: https://github.com/nextstrain/augur/pull/766
[#772]: https://github.com/nextstrain/augur/pull/772

## 13.0.0 (17 August 2021)

### Major Changes

* filter: Skip metadata records with ambiguous month information in the `date` column when grouping by month instead of randomly generating month values for those records. This change alters the behavior of the `filter` command for metadata with ambiguous month values. For these data, consider using `--group-by year` instead of `--group-by year month`. [#761][] (@huddlej)

### Features

* filter: When grouping by year or month, report the number of strains skipped due to ambiguous year and month both in the summary report at the end of filtering and in the `--output-log` contents [#761][] (@huddlej)

[#761]: https://github.com/nextstrain/augur/pull/761

## 12.1.1 (13 August 2021)

### Bug Fixes

* filter: Fix parsing of missing data in metadata [#758][] (@huddlej)
* filter: Fix probabilistic sampling with small values [#759][] (@huddlej)

[#758]: https://github.com/nextstrain/augur/pull/758
[#759]: https://github.com/nextstrain/augur/pull/759

## 12.1.0 (12 August 2021)

### Features

* export: Add support for custom legend and color scale specifications in Auspice config files [#727][] (@jameshadfield)
* utils: Add support for compressed strain name files (e.g., "include.txt.gz") [#730][] (@benjaminotter)
* filter: Rewrite internal logic to use pandas DataFrames ([#743][]), define filters and subsampling logic as individual functions ([#745][] and [#746][]), and iterate through chunks of metadata instead of loading all records into memory at once ([#750][]) (@tsibley, @huddlej)

### Bug Fixes

* distance: Change numeric type of distance output to float [#729][] (@benjaminotter)
* filter: Disable probabilistic sampling when users provide `--sequences-per-group` [#737][] (@benjaminotter)
* export: Provide correct missing file error messages for metadata and node data JSON inputs [#752][] (@benjaminotter)

[#727]: https://github.com/nextstrain/augur/pull/727
[#729]: https://github.com/nextstrain/augur/pull/729
[#730]: https://github.com/nextstrain/augur/pull/730
[#737]: https://github.com/nextstrain/augur/pull/737
[#743]: https://github.com/nextstrain/augur/pull/743
[#745]: https://github.com/nextstrain/augur/pull/745
[#746]: https://github.com/nextstrain/augur/pull/746
[#750]: https://github.com/nextstrain/augur/pull/750
[#752]: https://github.com/nextstrain/augur/pull/752

## 12.0.0 (13 April 2021)

### Major Changes

* filter: Date bounds (`--min-date` and `--max-date`) are now inclusive instead of exclusive such that records matching the given dates will pass date filters [#708][] (@benjaminotter)

### Bug Fixes

* refine: Recommend an alternate action when skyline optimization fails [#712][] (@huddlej)

### Features

* distance: Count insertion/deletion events once in pairwise distances [#698][] (@huddlej, @benjaminotter)
* distance: Optionally ignore specific list of characters defined in a distance map's top-level `ignored_characters` list [#707][] (@benjaminotter)
* filter: Allow `--subsample-max-sequences` without `--group-by` [#710][] (@benjaminotter)
* tree: Prefer `iqtree2` binary over `iqtree` when possible [#711][] (@benjaminotter)

[#698]: https://github.com/nextstrain/augur/pull/698
[#707]: https://github.com/nextstrain/augur/pull/707
[#708]: https://github.com/nextstrain/augur/pull/708
[#710]: https://github.com/nextstrain/augur/pull/710
[#711]: https://github.com/nextstrain/augur/pull/711
[#712]: https://github.com/nextstrain/augur/pull/712

## 11.3.0 (19 March 2021)

### Bug Fixes

* filter: Clarify how the `--priority` input affects subsampling in the command line help documentation [#695][]
* tests: Clean up outputs created by tests [#703][], ignore log files [#701][], and remove unnecessary Conda environment file [#702][]

### Features

* io: Add new `io` module with `open_file`, `read_sequences`, and `write_sequences` functions that support compressed inputs and outputs [#652][]
* parse, index, filter, mask: Add support for compressed inputs/outputs [#652][]
* export v2: Add optional `data_provenance` field to auspice JSON output for better provenance reporting in Auspice [#705][]

[#652]: https://github.com/nextstrain/augur/pull/652
[#695]: https://github.com/nextstrain/augur/pull/695
[#701]: https://github.com/nextstrain/augur/pull/701
[#702]: https://github.com/nextstrain/augur/pull/702
[#703]: https://github.com/nextstrain/augur/pull/703
[#705]: https://github.com/nextstrain/augur/pull/705

## 11.2.0 (8 March 2021)

### Bug Fixes

* ancestral: Mask positions that are ambiguous in all tip sequences before inferring ancestral sequence states, to avoid assigning arbitrary ancestral values based on rounding errors [#682][]
* titers: Add missing `kwargs` attribute to `TiterCollection` class [#690][]

### Documentation

* Update API documentation to include newer Python modules and the `index` subcommand [#687][]
* Remove Zika and TB tutorials in favor of copies in docs.nextstrain.org [#689][]

### Features

* filter: Enable filtering by metadata only such that sequence inputs/outputs are optional and metadata/strain list outputs are now possible [#679][]
* filter: Enable extraction of sequences from multiple lists of strains with a new `--exclude-all` flag and support for multiple inputs to the `--include` argument [#679][]

[#679]: https://github.com/nextstrain/augur/pull/679
[#682]: https://github.com/nextstrain/augur/pull/682
[#687]: https://github.com/nextstrain/augur/pull/687
[#689]: https://github.com/nextstrain/augur/pull/689
[#690]: https://github.com/nextstrain/augur/pull/690

## 11.1.2 (16 February 2021)

### Bug Fixes

* index: Remove call to deprecated BioPython SeqIO.close method [#684][]

[#684]: https://github.com/nextstrain/augur/pull/684

## 11.1.1 (16 February 2021)

### Bug Fixes

* filter: Retry probabilistic subsampling when it doesn't select any samples [#676][]
* titers: Skip infinite log titer values caused by zero-valued raw titers [#677][]
* setup: Include license file with distribution artifacts instead of Python installation root [#678][]

[#676]: https://github.com/nextstrain/augur/pull/676
[#677]: https://github.com/nextstrain/augur/pull/677
[#678]: https://github.com/nextstrain/augur/pull/678

## 11.1.0 (12 February 2021)

### Bug Fixes

* align/utils: Improve external command errors [#638][]
* filter: Fix parsing of priorities files to allow spaces in sequence IDs [#668][]
* utils: Ensure columns `strain` and `name` in metadata get parsed as strings [#669][]

### Features

* index/filter: Add new `index` subcommand and optional `--sequence-index` argument for `filter` command to enable filtering without inspecting sequences [#651][]
* titers: Bump supported cvxopt version to latest 1.* release [#671][]

[#638]: https://github.com/nextstrain/augur/pull/638
[#651]: https://github.com/nextstrain/augur/pull/651
[#668]: https://github.com/nextstrain/augur/pull/668
[#669]: https://github.com/nextstrain/augur/pull/669
[#671]: https://github.com/nextstrain/augur/pull/671

## 11.0.0 (22 January 2021)

### Major Changes

* filter: Use probabilistic sampling by default when requesting a maximum number of sequences to subsample with `--subsample-max-sequences`. Adds `--no-probabilistic-sampling` flag to disable this default behavior and prevent users from requesting fewer maximum sequences than there are subsampling groups. [#659][]

[#659]: https://github.com/nextstrain/augur/pull/659

## 10.3.0 (14 January 2021)

### Bug Fixes

* scripts: Fix typo in `verify_meta_json.py` [#656][] (@felixonmars)
* CI: Use the expected Python version in conda environments [#658][]
* CI: Minimize codecov feedback [#661][]

### Features

* frequencies: Add `--pivot-interval-units` argument with support for weekly pivots [#660][]
* frequencies: Add support for ISO dates for `--min-date` and `--max-date` arguments [#660][]

[#656]: https://github.com/nextstrain/augur/pull/656
[#658]: https://github.com/nextstrain/augur/pull/658
[#660]: https://github.com/nextstrain/augur/pull/660
[#661]: https://github.com/nextstrain/augur/pull/661

## 10.2.0 (1 January 2021)

### Features

* filter: Add `--probablistic-sampling` flag to allow subsampling with `--subsample-max-sequences` when the number of groups exceeds the requested number of samples [#629][]
* scripts: Add script to identify emerging clades from existing Nextstrain build JSONs [#653][]
* docs: Add instructions to update conda installations prior to installing Augur [#655][]

[#629]: https://github.com/nextstrain/augur/pull/629
[#653]: https://github.com/nextstrain/augur/pull/653
[#655]: https://github.com/nextstrain/augur/pull/655

## 10.1.1 (16 November 2020)

### Bug Fixes

* dependencies: Require the most recent minor versions of TreeTime (0.8.X) to fix numpy matrix errors [#633][]

[#633]: https://github.com/nextstrain/augur/pull/633

## 10.1.0 (13 November 2020)

### Features

* docs: Migrate non-reference documentation to docs.nextstrain.org [#620][]
* filter: Add `--exclude-ambiguous-dates-by` flag to enable exclusion of samples with ambiguous dates [#623][] and [#631][]

[#620]: https://github.com/nextstrain/augur/pull/620
[#623]: https://github.com/nextstrain/augur/pull/623
[#631]: https://github.com/nextstrain/augur/pull/631

## 10.0.4 (6 November 2020)

### Bug Fixes

* tree: Use a more generic approach to escape special characters from alignment sequence names prior to running IQ-TREE [#625][]
* filter: Reduce memory usage by not reading sequences into memory [#627][]

[#625]: https://github.com/nextstrain/augur/pull/625
[#627]: https://github.com/nextstrain/augur/pull/627

## 10.0.3 (23 October 2020)

### Bug Fixes

* refine: Report divergence by number of mutations as an integer instead of a floating point value [#618][]
* validate: Allow internal nodes with a single child and do not allow duplicate node names [#621][]

[#618]: https://github.com/nextstrain/augur/pull/618
[#621]: https://github.com/nextstrain/augur/pull/621

## 10.0.2 (8 September 2020)

### Bug Fixes

* align: Remove references to BioPython's deprecated `Alphabet` attributes [#615][]
* Pin BioPython dependency to a max supported version to prevent breaking changes to augur in the future [#615][]

[#615]: https://github.com/nextstrain/augur/pull/615

## 10.0.1 (8 September 2020)

### Bug Fixes

* ancestral: Clarify default values for inference of ambiguous bases [#613][]

[#613]: https://github.com/nextstrain/augur/pull/613

## 10.0.0 (17 August 2020)

### Major Changes

* Remove Snakemake as a dependency of the augur Python package [#557][]
  * Updated documentation to reflect external Snakemake dependency [#600][] and Snakemake's required `--cores` argument [#599][]
* utils: `read_colors` refactor [#588][]
  * raises an exception when the requested color file is missing instead of printing a warning to stdout
  * splits out logic to parse colors file into separate classes (`util_support/color_parser.py` and `util_support/color_parser_line.py`) with unit tests
* utils: `read_metadata` interface improvements
  * raises exceptions when 1) input file is missing or unreadable or 2) required columns (`strain` or `name`) are missing instead of failing silently [#584][]
  * automatically detects delimiter in metadata file instead of assuming delimiter based on filename extension [#587][]
* utils: `read_node_data` interface improvements [#595][], [#605][]
  * exits with a nonzero code when node data node names don't match tree nodes and when the input tree cannot be loaded
  * refactors logic to read node data into separate classes with unit tests

### Bug Fixes

* ancestral: Fix docstring for `collect_mutations_and_sequences` [4c474a9][]
* parse: Fix date parsing bug caused by a change in the API for `parse_time_string` in pandas 1.1.0 [#601][]
* refine: Enable divergence unit scaling without timetree [e9b3eec][]
* tree: Use IQ-TREE's `-nt AUTO` mode when users request more threads than there are input sequences, avoiding an IQ-TREE error [#598][]

### Features

* docs: Document support for installation from Bioconda [#604][]
* filter: Add `--subsample-max-sequences` argument to limit the maximum number of sequences to be included in subsampled output [#593][]
* mask: Add `--mask-invalid` flag to mask invalid nucleotides from FASTA files [#592][]

[#557]: https://github.com/nextstrain/augur/pull/557
[#584]: https://github.com/nextstrain/augur/pull/584
[#587]: https://github.com/nextstrain/augur/pull/587
[#588]: https://github.com/nextstrain/augur/pull/588
[#592]: https://github.com/nextstrain/augur/pull/592
[#593]: https://github.com/nextstrain/augur/pull/593
[#595]: https://github.com/nextstrain/augur/pull/595
[#598]: https://github.com/nextstrain/augur/pull/598
[#599]: https://github.com/nextstrain/augur/pull/599
[#600]: https://github.com/nextstrain/augur/pull/600
[#601]: https://github.com/nextstrain/augur/pull/601
[#604]: https://github.com/nextstrain/augur/pull/604
[#605]: https://github.com/nextstrain/augur/pull/605
[e9b3eec]: https://github.com/nextstrain/augur/commit/e9b3eec670b9603874e195cc1ccd4f3c1aeef5dd
[4c474a9]: https://github.com/nextstrain/augur/commit/4c474a96232e9cc333e3fc4c0971336a090b703c

## 9.0.0 (29 June 2020)

### Major Changes

* align: The API to the `read_sequences` function now returns a list of sequences instead of a dictionary [#536][]

### Bug Fixes

* align: Prevent duplicate strains warning when using `--reference-name` [#536][]
* docs: Sync and deduplicate installation documentation from README to main docs [#578][]
* export: Flexibly disambiguate multiple publications by the same author [#581][]
* frequencies: Avoid interpolation of a single data point during frequency estimation with sparse data [#569][]
* parse: Actually remove commas during prettify when this behavior is requested [#573][]
* tests: Always use the local helper script (`bin/augur`) to run tests instead of any globally installed augur executables [#527][]
* tree: Keep log files after trees are built [#572][]
* utils: Do not attempt to parse dates with only ambiguous months (e.g., 2020-XX-01) [#532][]
* utils: Parse `name` column of metadata as a data field instead of a pandas DataFrame attribute [#564][]

### Features

* docs: Updates description of how missing data are handled by `augur traits`
* filter: Add support for ISO 8601 dates (YYYY-MM-DD) for `--min-date` and `--max-date` [#568][]
* tests: Add tests for utilities (ambiguous date parsing [#532][] and `run_shell_command` [#577][]), parse [#573][], and translate [#546][]
* tree: Allow VCF input without an `--exclude-sites` argument [#565][]

[#527]: https://github.com/nextstrain/augur/pull/527
[#532]: https://github.com/nextstrain/augur/pull/532
[#536]: https://github.com/nextstrain/augur/pull/536
[#546]: https://github.com/nextstrain/augur/pull/546
[#564]: https://github.com/nextstrain/augur/pull/564
[#565]: https://github.com/nextstrain/augur/pull/565
[#568]: https://github.com/nextstrain/augur/pull/568
[#569]: https://github.com/nextstrain/augur/pull/569
[#572]: https://github.com/nextstrain/augur/pull/572
[#573]: https://github.com/nextstrain/augur/pull/573
[#577]: https://github.com/nextstrain/augur/pull/577
[#578]: https://github.com/nextstrain/augur/pull/578
[#581]: https://github.com/nextstrain/augur/pull/581

## 8.0.0 (8 June 2020)

### Major Changes

* utils: Add a consolidated generic `load_mask_sites` function and specific `read_mask_file` and `read_bed_file` functions for reading masking sites from files. Changes the Python API by moving mask-loading functionality out of augur mask and tree into utils [#514][] and [#550][]
* mask: Parse BED files as zero-indexed, half-open intervals [#512][]

### Bug Fixes

* traits: Export mugration models to the same output directory as traits JSON [#544][]
* Explicitly open files with UTF-8 file encoding [#499][], [#503][], and [#560][]
* refine: Only request confidence intervals from TreeTime when no clock rate is provided [#548][]
* refine: Catch failed skyline optimization [#556][]

### Features

* align: Report insertions stripped during alignment [#449][]
* Require minimum pandas version of 1.0.0 [#488][]
* parse: Reduce memory use and clarify code with standard Python idioms [#496][]
* mask: Allow masking of specific sites passed by the user with `--mask-sites` and masking of a fixed number of sites from the beginning or end of each sequence with `--mask-from-beginning` and `--mask-from-end` [#512][]
* clades, import: Use `defaultdict` to simplify code [#533][]
* tests: Add initial functional tests of the augur command line interface using Cram [#542][]
* refine: Add a `--seed` argument to set the random seed for more reproducible outputs across runs [#542][]
* ancestral, refine, and traits: Print the version of TreeTime being used for these commands [#552][]
* filter: Add support for flexible pandas-style queries with new `--query` argument [#555][]
* export: Allow display defaults for transmission lines [#561][]

[#449]: https://github.com/nextstrain/augur/pull/449
[#488]: https://github.com/nextstrain/augur/pull/488
[#496]: https://github.com/nextstrain/augur/pull/496
[#499]: https://github.com/nextstrain/augur/pull/499
[#503]: https://github.com/nextstrain/augur/pull/503
[#512]: https://github.com/nextstrain/augur/pull/512
[#514]: https://github.com/nextstrain/augur/pull/514
[#533]: https://github.com/nextstrain/augur/pull/533
[#542]: https://github.com/nextstrain/augur/pull/542
[#544]: https://github.com/nextstrain/augur/pull/544
[#548]: https://github.com/nextstrain/augur/pull/548
[#550]: https://github.com/nextstrain/augur/pull/550
[#552]: https://github.com/nextstrain/augur/pull/552
[#555]: https://github.com/nextstrain/augur/pull/555
[#556]: https://github.com/nextstrain/augur/pull/556
[#560]: https://github.com/nextstrain/augur/pull/560
[#561]: https://github.com/nextstrain/augur/pull/561

## 7.0.2 (7 April 2020)

### Bug Fixes

* filter: Fix regression introduced in 7.0.0 which caused an error to be raised
  if a priorities file didn't include every sequence.  Sequences which are not
  explicitly listed will once again default to a priority of 0. [#530][]

[#530]: https://github.com/nextstrain/augur/pull/530


## 7.0.1 (7 April 2020)

### Bug Fixes

* Fix typo with Python classifiers in setup.py

## 7.0.0 (7 April 2020)

### Major Changes

* Drop support for Python 3.4 and 3.5 [#482][]
* Drop support for `--output` flag in augur ancestral, clades, sequence-traits, traits, and translate in favor of `--output-node-data` flag [#529][]

### Features

* improve testing by
	* adding a simple shell script to run tests and improving pytest configuration and output [#463][]
	* adding code coverage reports ([#486][], [#491][]) and integration with codecov.io [#508][]
	* adding unit tests for align ([#477][]), filter ([#478][], [#479][], [#487][]), utils ([#501][])
* align: reverse complement sequences when necessary using mafft’s autodirection flag [#467][]
* align: speed up replacement of gaps with “ambiguous” bases [#474][]
* mask: add support for FASTA input files [#493][]
* traits: bump TreeTime version to 0.7.4 and increase maximum number of unique traits allowed from 180 to 300 [#495][]

### Bug Fixes

* align: enable filling gaps in input sequences even if no reference is provided instead of throwing an exception [#466][]
* align: detect duplicate sequences by comparing sequence objects instead of (often truncated) string representations of those objects [#468][]
* import_beast: use raw strings for regular expressions to avoid syntax errors in future versions of Python [#469][]
* scripts: update exception syntax to new style [#484][]
* filter: fail loudly when a given priority file is invalid and exit instead of just printing an error [#487][]

### Documentation

* README: document need for brew tap prior to brew install [#471][]
* DEV_DOCS: add a proper section on how to write and run tests [#481][]
* README: improve discoverability of core documentation with RTD badge [#490][]

[#463]: https://github.com/nextstrain/augur/pull/463/
[#466]: https://github.com/nextstrain/augur/pull/466
[#467]: https://github.com/nextstrain/augur/pull/467
[#468]: https://github.com/nextstrain/augur/pull/468
[#469]: https://github.com/nextstrain/augur/pull/469
[#471]: https://github.com/nextstrain/augur/pull/471
[#474]: https://github.com/nextstrain/augur/pull/474
[#477]: https://github.com/nextstrain/augur/pull/477
[#478]: https://github.com/nextstrain/augur/pull/478
[#479]: https://github.com/nextstrain/augur/pull/479
[#481]: https://github.com/nextstrain/augur/pull/481
[#482]: https://github.com/nextstrain/augur/pull/482
[#484]: https://github.com/nextstrain/augur/pull/484
[#486]: https://github.com/nextstrain/augur/pull/486
[#487]: https://github.com/nextstrain/augur/pull/487
[#490]: https://github.com/nextstrain/augur/pull/490
[#491]: https://github.com/nextstrain/augur/pull/491
[#493]: https://github.com/nextstrain/augur/pull/493
[#495]: https://github.com/nextstrain/augur/pull/495
[#501]: https://github.com/nextstrain/augur/pull/501
[#508]: https://github.com/nextstrain/augur/pull/508
[#529]: https://github.com/nextstrain/augur/pull/529

## 6.4.3 (25 March 2020)

### Bug Fixes

* align: Remove reference sequence from alignments even when no gaps exist in
  input sequences relative to the reference. Thank you @danielsoneg!
  [#456][]

### Documentation

* Reorganize README, improve findability of documentation, and add separate dev
  docs.
  [#461][]

[#456]: https://github.com/nextstrain/augur/pull/456
[#461]: https://github.com/nextstrain/augur/pull/461

## 6.4.2 (17 March 2020)

### Bug Fixes

* Require Snakemake less than 5.11 to avoid a breaking change.  The `--cores`
  argument is now required by 5.11, which will affect many existing augur-based
  workflows.  Reported upstream as
  [snakemake/snakemake#283](https://github.com/snakemake/snakemake/issues/283).

* align: Run mafft with the `--nomemsave` option.  This makes alignments of
  sequences over 10k in length run much, much faster in the general case and
  shouldn't cause issues for most modern hardware.  We may end up needing to
  add an off-switch for this mode if it causes issues for other users of augur,
  but the hope is that it will make things just magically run faster for most
  folks!  There is likely more tuning that could be done with mafft, but this
  is a huge improvement in our testing.
  [#458][]

* align: Ignore blank lines in `--include` files.  Thanks @CameronDevine!
  [#451][]

* align: Properly quote filenames when invoking mafft.  Thanks @CameronDevine!
  [#452][]


[#451]: https://github.com/nextstrain/augur/pull/451
[#452]: https://github.com/nextstrain/augur/pull/452
[#458]: https://github.com/nextstrain/augur/pull/458


## 6.4.1 (4 March 2020)

### Bug Fixes

* export: AA labels are now exported for branches where a clade is also labeled
  [See PR 447](https://github.com/nextstrain/augur/pull/447)

* export / validation: a dataset title is no longer required

* release script now works on MacOS & code-signing is optional
  [See PR 448](https://github.com/nextstrain/augur/pull/448)

* traits: Missing data is correctly handled

## 6.4.0 (26 February 2020)

### Features

* align: New sequences can now be added to an existing alignment.  [#422][]

* align: Multiple sequence files can be provided as input. [#422][]

* align: Extra debugging files such as `*.pre_aligner.fasta` and
  `*.post_aligner.fasta` are no longer produced by default.  To request them,
  pass the `--debug` flag. [#422][]

* align: De-duplicate input sequences, with a warning. [#422][]

* export v2: Add support for the `branch_label` property in `display_defaults`,
  which was recently added to Auspice. [#445][]

### Bug fixes

* align: Exits with an error earlier if arguments are invalid instead of only
  printing a warning. [#422][]

* align: Performs more error checking and clarifies the help and error
  messages. [#422][]

* export v2: Traits which are filters but not colorings are now exported as
  well, instead of being left out. [#442][]

* export v2: Exits non-zero when validation fails, instead of masking errors.
  [#441][]

* validate: In order to improve clarity, messages now include the filenames
  involved and distinguish between schema validation and internal consistency
  checks. [#441][]

[#422]: https://github.com/nextstrain/augur/pull/422
[#441]: https://github.com/nextstrain/augur/pull/441
[#442]: https://github.com/nextstrain/augur/pull/442
[#445]: https://github.com/nextstrain/augur/pull/445


## 6.3.0 (13 February 2020)

### Features

* Augur `refine`, `ancestral` and `traits` now use the
  [upgraded TreeTime v0.7](https://github.com/neherlab/treetime/releases/tag/v0.7.0)
  This should have a number of under-the-hood improvements.
  [See PR 431](https://github.com/nextstrain/augur/pull/431)
* ancestral: New options to either `--keep-ambiguous` or `--infer-ambiguous`. If using
  `--infer-ambiguous` the previous behavior will be maintained in which tips with `N` will have
  their nucleotide state inferred. If using `--keep-ambiguous`, these tips will be left as `N`.
  With this upgrade, we are still defaulting to `--infer-ambiguous`, however, we plan to swap
  default to `--keep-ambiguous` in the future. If this distintion matters to you, we would suggest
  that you explicitly record `--keep-ambiguous` / `--infer-ambiguous` in your build process.
  [Also part of PR 431](https://github.com/nextstrain/augur/pull/431)
* traits: Allow input of `--weights` which references a `.tsv` file in the following format:
  ```
  division	Hubei	10.0
  division	Jiangxi	1.0
  division	Chongqing	1.0
  ```
  where these weights represent equilibrium frequencies in the CTMC transition model. We imagine the
  primary use of user-specified weights to correct for strong sampling biases in available data.
  [See PR 443](https://github.com/nextstrain/augur/pull/443)

### Bug fixes

* Improvements to make shell scripts run more easily on Windows.
  [See PR 437](https://github.com/nextstrain/augur/pull/437)

## 6.2.0 (25 January 2020)

### Features

* refine: Include `--divergence-units` option to distinguish between `mutations`
  and `mutations-per-site`. Keep `mutations-per-site` as default behavior.
  [See PR 435](https://github.com/nextstrain/augur/pull/435)

### Bug fixes

* utils: Support v2 auspice JSONs in json_to_tree utility function.
  [See PR 432](https://github.com/nextstrain/augur/pull/432)

## 6.1.1 (17 December 2019)

### Bug fixes

* frequencies: Fix bug in string matching for weighted frequencies introduced in
  v6.1.0. [See PR 426](https://github.com/nextstrain/augur/pull/426).

## 6.1.0 (13 December 2019)

### Features

* export: Include `--description` option to pass in a Markdown file with dataset
  description. This is displays in Auspice in the footer. For rationale,
  [see Auspice issue 707](https://github.com/nextstrain/auspice/issues/707) and
  for Augur changes [see PR 423](https://github.com/nextstrain/augur/pull/423).

### Bug fixes

* frequencies: Fix weighted frequencies when weight keys are unrepresented.
  [See PR 420](https://github.com/nextstrain/augur/pull/420).


## 6.0.0 (10 December 2019)

### Overview

Version 6 is a major release of augur affecting many augur commands. The format
of the exported JSON (v2) has changed and now merges the previously separate
files containing tree and meta information. To maintain backward compatibility,
the export command was split into `export v1` (old) and `export v2` (new).
Detailed release notes are provided in the augur documentation [on
read-the-docs](https://nextstrain-augur.readthedocs.io/en/stable/releases/v6.html).
For a migration guide, consult
[migrating-v5-v6](https://nextstrain-augur.readthedocs.io/en/stable/releases/migrating-v5-v6.html).

### Major features / changes

* export: Swap from a separate `_tree.json` and `_meta.json` to a single
  "unified" `dataset.json` output file
* export: Include additional command line options to alleviate need for Auspice
  config
* export: Include option for reference sequence output
* export: Move to GFF-style annotations
* export: Validate exported JSONs against schema
* ancestral: Allow output of FASTA and JSON files
* import: Include `import beast` command to import labeled BEAST MCC tree
* parse: Include `--prettify-fields` option to cleanup metadata fields
* Documentation improvements

### Minor features / changes

* colors.tsv: Allow whitespace, but insist on tab delimiting
* lat_longs.tsv: Allow whitespace, but insist on tab delimiting
* Remove code for old "non-modular" augur, old "non-modular" builds and Python
  tests
* Improve test builds
* filter: More interpretable output of how many sequences have been filtered
* filter: Additional flag `--subsample-seed` to seed the random number generator
  and thereby make subsampling reproducible
* sequence-traits: Numerical output as originally intended, but required an
  Auspice bugfix
* traits: Explanation of what is considered missing data & how it is interpreted
* traits: GTR models are exported in the output JSON for better accountability &
  reproducibility


## 5.4.1 (12 November 2019)

### Bug fixes

* export v1: Include `--minify-json` option that was mistakenly not included in PR 398.
  [See PR 409](https://github.com/nextstrain/augur/pull/409)

## 5.4.0 (7 November 2019)

### Features

* frequencies: Include `--minimal-clade-size-to-estimate` command line option.
  [See PR 383](https://github.com/nextstrain/augur/pull/383)
* lbi: Include `--no-normalization` command line option.
  [See PR 380](https://github.com/nextstrain/augur/pull/380)

### Compatibility fixes

* export: Include `v1` subcommand to allow forwards compatibiliy with Augur v6 builds.
  [See PR 398](https://github.com/nextstrain/augur/pull/398)

### Bug fixes

* export: Include warning if using a mismatched v6 translate file.
  [See PR 392](https://github.com/nextstrain/augur/pull/392)
* frequencies: Fix determination of interval for clipping of non-informative pivots

## 5.3.0 (9 September 2019)

### Features

* export: Improve printing of error messages with missing or conflicting author
  data. [See issue 274](https://github.com/nextstrain/augur/issues/274)
* filter: Improve printing of dropped strains to include reasons why strains were
  dropped. [See PR 367](https://github.com/nextstrain/augur/pull/367)
* refine: Add support for command line flag `--keep-polytomies` to not resolve
  polytomies when producing a time tree.
  [See PR 345](https://github.com/nextstrain/augur/pull/345)

### Bug fixes

* Catch and throw error when there are duplicate strain names.
  [See PR 356](https://github.com/nextstrain/augur/pull/356)
* Fix missing annotation of "parent" attribute for the root node
* Run shell commands with more robust error checking.
  [See PR 350](https://github.com/nextstrain/augur/pull/350)
* Better handling of rerooting options for trees without temporal information.
  [See issue 348](https://github.com/nextstrain/augur/issues/348)

### Data

* Small fixes in geographic coordinate file

## 5.2.1 (4 August 2019)

### Bug fixes

* Print more useful error message if Python recursion limit is reached.
  [See issue 328](https://github.com/nextstrain/augur/issues/328)
* Print more useful error message if vcftools if missing.
  [See PR 312](https://github.com/nextstrain/augur/pull/321)

### Development

* Significantly relax version requirements specified in setup.py for biopython,
  pandas, etc... Additionally, move lesser used packages (cvxopt, matplotlib,
  seaborn) into an "extras_require" field. This should reduce conflicts with
  other pip installed packages.
  [See PR 323](https://github.com/nextstrain/augur/pull/323)

### Data

* Include additional country lat/longs in base data

## 5.2.0 (23 July 2019)

### Features

* ancestral: Adds a new flag `--output-sequences` and logic to support saving
  ancestral sequences and leaves from the given tree to a FASTA file. Also adds a
  redundant, more specific flag `--output-node-data` that will replace the current
  `--output` flag in the next major version release of augur. For now, we issue a
  deprecation warning when the `--output` flag is used. Note that FASTA output is
  only allowed for FASTA inputs and not for VCFs. We don't allow FASTA output for
  VCFs anywhere else and, if we did here, the output files would be very large.
  [See PR 293](https://github.com/nextstrain/augur/pull/293)

* frequencies: Allow `--method kde` flag to compute frequencies via KDE kernels.
  This complements existing method of `--method diffusion`. Generally, KDE
  frequencies should be more robust and faster to run, but will not project as
  well when forecasting frequencies into the future.
  [See PR 271](https://github.com/nextstrain/augur/pull/271)

### Bug fixes

* ancestral, traits, translate: Print warning if supplied tree is missing internal
  node names (normally provided by running `augur refine`).
  [See PR 283](https://github.com/nextstrain/augur/pull/283)

* Include pip in Conda enviroment file.
  [See PR 309](https://github.com/nextstrain/augur/pull/309)

### Documentation

* Document environment variables respected by Augur

### Development

* Remove matplotlib and seaborn from `setup.py` install. These are still called a
  few places in augur (like `titers.validate()`), but it was deemed rare enough
  that remove this from `setup.py` would ease general install for most users.
  Additionally, the ipdb debugger has been moved to dev dependencies.
  [See PR 291](https://github.com/nextstrain/augur/pull/291)

* Refactor logic to read trees from multiple formats into a function. Adds a new
  function `read_tree` to the `utils` module that tries to safely handle reading
  trees in multiple input formats.
  [See PR 310](https://github.com/nextstrain/augur/pull/310)


## 5.1.1 (1 July 2019)

### Features

* tree: Add support for the GTR+R10 substitution model.
* tree: Support parentheses in node names when using IQ-TREE.

### Bug fixes

* Use the center of the UK for its coordinates instead of London.
* filter: Mark `--output` required, which it always was but wasn't marked.
* filter: Avoid error when no excluded strains file is provided.
* export: Fix for preliminary version 2 schema support.
* refine: Correct error handling when the tree file is missing or empty.

### Documentation

* Add examples of Augur usage in the wild.
* Rename and reorganize CLI and Python API pages a little bit to make "where do
  I start learning to use Augur?" clearer to non-devs.

### Development

* Relax version requirements of pandas and seaborn.  The hope is this will make
  installation smoother (particularly alongside other packages which require
  newer pandas versions) while not encountering breaking changes in newer
  versions ourselves.


## 5.1.0 (29 May 2019)

### Documentation

* Documentation is now available online for the augur CLI and Python API via
  Read The Docs: <https://nextstrain-augur.readthedocs.io>.  The _latest_
  version on RTD points to the git master branch, and the _stable_ version to
  the most recent tagged release.  Instructions for building the docs locally
  are in the README.


## 5.0.0 (26 May 2019)

### Features

* ancestral: New option to `--keep-ambiguous`, which will not infer nucleotides at
  ambiguous (N) sites on tip sequences and instead leave as 'N'
  [See PR 280.](https://github.com/nextstrain/augur/pull/280)
* ancestral: New option to `--keep-overhangs`, which will not infer nucleotides for
  gaps on either side of the alignment and instead leave as '-'.
  [See PR 286.](https://github.com/nextstrain/augur/pull/286)
* clades: This module has been reconfigured to identify clade defining mutations on
  top of a reference rather than identifying mutations along the tree. The command
  line arguments are the same except for the addition of `--reference`, which
  explicitly passes in a reference sequence. If `--reference` is not defined, then
  reference will be drawn from the root node of the phylogeny by looking for
  `sequence` attribute attached to root node of `--tree`.
  [See PR 288.](https://github.com/nextstrain/augur/pull/288)
* refine: Revise rooting behavior. Previously `--root` took 'best', 'residual', 'rsq'
  and 'min_dev' as options. In this update `--root` takes 'best', least-squares',
  'min_dev' and 'oldest' as rooting options. This eliminates 'residual' and 'rsq'
  as options. This is a **backwards-incompatible** change. This requires updating
  TreeTime to version 0.5.4 or above.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* refine: Add `--keep-root` option that overrides `--root` specification to preserve
  tree rooting.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* refine: Add `--covariance` and `--no-covariance` options that specify TreeTime
  behavior.
  [See PR 263.](https://github.com/nextstrain/augur/pull/263)
* titers: This command now throws an `InsufficientDataException` if there are not
  sufficient titers to infer a model. This is paired with a new `--allow-empty-model`
  flag that proceeds past the `InsufficientDataException` and writes out a model
  JSON corresponding to an 'empty' model.
  [See PR 281.](https://github.com/nextstrain/augur/pull/281)
* By default JSONs are written with `index=1` to give a pretty-printed JSON. However,
  this adds significant file size to large tree JSONs. If the environment variable
  `AUGUR_MINIFY_JSON` is set then minified JSONs are printed instead. This mirror the
  explicit `--minify-json` argument available to `augur export`.
  [See PR 278.](https://github.com/nextstrain/augur/pull/278)

### Bug fixes

* export: Cast numeric values to strings for export.
  [See issue 287.](https://github.com/nextstrain/augur/issues/287)
* export: Legend order preserves ordering passed in by user for traits that have
  default colorings ('country' and 'region').
  [See PR 284.](https://github.com/nextstrain/augur/pull/284)
* refine: Previously, the `--root` argument was silently ignored when no timetree
  was inferred. Re-rooting with an outgroup is sensible even without a timetree.
  [See PR 282.](https://github.com/nextstrain/augur/pull/282)

## 4.0.0 (24 April 2019)

### Features

* distance: New interface for specifying distances between sequences. This is
  a **backwards-incompatible** change. Refer to `augur distance --help` for
  all the details.

* export: Add a `--minify-json` flag to omit indentation in Auspice JSONs.

### Bug fixes

* frequencies: Emit one-based coordinates (instead of zero-based) for KDE-based
  mutation frequencies

### Data

* Include additional country lat/longs in base data


## 3.1.8 (13 February 2019)

### Bug fixes

* titers: fix calculation of `mean_potentency` for model export

## 3.1.7 (5 February 2019)

### Bug fixes

* Update to TreeTime 0.5.3
* tree: Fix bug in printing causing errors in Python versions <3.6
* tree: Alter site masking to not be so memory intensive

## 3.1.6 (29 January 2019)

### Features

* filter: Allow negative matches to `--exclude-where`. For example,
 `--exclude-where country!=usa` would exclude all samples where metadata `country` does
 not equal `usa`.
* tree: Allow `--exclude-sites` to work with FASTA input. Ensure that indexing of input
 sites is one-based.

### Bug fixes

* fix loading of strains when loading titers from file, previously strains had not been
 filtered to match the tree appropriately

## 3.1.5 (13 January 2019)

### Features

* frequencies: Add `--ignore-char` and `--minimal-clade-size` as options.
* frequencies: Include `--stiffness` and `--inertia` as options.
* titers: Allow multiple titer date files in `--titers` import.

### Bug fixes

* filter: Fix `--non-nucleotide` call to include `?` as allowed character.
* tree: Fix `--method raxml` to properly delimit interim RAxML output so that
  simultaneous builds don't conflict.

### Data

* Include additional country lat/longs in base data

## 3.1.4 (1 January 2019)

### Bug fixes

* frequencies: Include `counts` in `augur frequencies` output JSON to support
  downstream plotting.

### Data

* Include additional country lat/longs in base data

## 3.1.3 (29 December 2018)

### Features

* filter: Add `--non-nucleotide` option to remove sequences with non-conforming
  nucleotide characters.

### Bug fixes

* Revise treatment of `-`, ` ` in `augur parse` to leave `-` as is and remove white
  space. Also delimit `[` and `]` to `_`.
* Fix bug in naming of temp IQTREE fixes to prevent conflicts from simultaneous builds.

### Data

* Include additional country lat/longs in base data

### Development

* Remove non-modular measles build in favor of [nextstrain/measles](https://github.com/nextstrain/measles) repo.

## 3.1.2 (21 December 2018)

### Bug fixes

* Update dependencies

## 3.1.1 (21 December 2018)

### Bug fixes

* filter: Fix `--include-where`. Adds an `all_seq` variable needed by the logic to
  include records by value. This was previously working for VCF but threw an exception
  for sequences in FASTA format.
* Update flu reference viruses and lat longs.
* Update dependencies

## 3.1.0 (18 December 2018)

### Features

* reconstruct-sequences: Include `augur reconstruct-sequences` module that reconstructs
  alignments from mutations inferred on the tree
* distance: Include `augur distance` module that calculates the distance between amino
  acid sequences across entire genes or at a predefined subset of sites
* lbi: Include `augur lbi` module that calculates local branching index (LBI) for a
  given tree and one or more sets of parameters.
* frequencies: Include `--method kde` as option to `augur frequencies`, separate from the
  existing `--method diffusion` logic. KDE frequencies are faster and better for smaller
  clades but don't extrapolate as well as diffusion frequencies.
* titers: Enable annotation of nodes in a tree from the substitution model

## 3.0.5.dev1 (26 November 2018)

### Bug fixes

* translate: Nucleotide ("nuc") annotation for non-bacterial builds starts at 0
  again, not 1, fixing a regression.

### Documentation

* Schemas: Correct coordinate system description for genome start/end
  annotations.


## 3.0.4.dev1 (26 November 2018)

### Bug fixes

* validate: Fix regression for gene names containing an asterisk.

### Development

* Fix Travis CI tests which were silently not running.


## 3.0.3.dev1 (26 November 2018)

### Features

* refine: Add a `--clock-std-dev` option

* traits: Add a `--sampling-bias-correction` option for mugration model

* validate: Gene names in tree annotations may now contain hyphens.  Compatible
  with Auspice version 1.33.0 and later.

* All JSON is now emitted with sorted keys, making it easier to diff and run
  other textual comparisons against output.

### Bug fixes

* filter: Only consider A, T, C, and G when calculating sequence length for the
  `--min-length` option.

* filter: Allow comments in files passed to `--exclude`.

* filter: Ignore case when matching trait values against excluded values.

* Normalize custom geographic names to lower case for consistent matching.

### Data

* Fix typo in geographic entry for `netherlands`.

* Schemas: Reconcile naming patterns used in gene definitions and tree
  annotations.

### Development

* Upgrade TreeTime dependency to 0.5.x and at least 0.5.1.

* Add an `environment.yml` file for use with `conda env create`.

* Stop testing under Python 2.7 on Travis CI.


## 3.0.2.dev1 (27 September 2018)

### Bug fixes

* translate: Fix broken `--help` message


## 3.0.1.dev1 (27 September 2018)

### Features

* align and tree: The --nthreads option now accepts the special value "auto" to
  automatically set the number of threads to the number of CPU cores available.

* Alias `augur --version` to `augur version`

### Bug fixes

* tree: The --nthreads option is now respected.  Previously all tree builders
  were ignoring the value and using either 2 threads (RAxML, IQ-TREE) or as
  many threads as cores (FastTree, if the OpenMP version).

* translate: Check for and, if necessary pad, nucleotide sequences which aren't
  a multiple of 3 earlier to avoid errors later.

* export: Optionally write inferred nucleotide and amino acid sequences (or
  mutations) to a separate file.

* export: Omit genes with no amino acid mutations.

* validate: Allow underscores in gene names.

* refine: Remove unused --nthreads argument.

* ancestral, filter, tree, refine: Exit 1 instead of -1 on error.

* Print the help message, instead of throwing an exception, when `augur` is run
  without arguments.

### Documentation

* Briefly describe each command in its `--help` output and in the global `augur
  --help` output.

* Revamp README to emphasize new, modular augur and make it suitable for
  inclusion on PyPi.

* Reconciled conflicting license declarations; augur is AGPLv3 (not MIT)
  licensed like the rest of Nextstrain.

* Include URLs for bug reports, the change log, and the source on PyPi.

### Data

* Geographic coordinates added for the Netherlands and the Philippines.

### Development

* Reset the `release` branch when rewinding a failed local release process.

* Refactor the augur program and command architecture for improved
  maintainability.


## 3.0.0.dev3 (4 September 2018)

### Development

* Use an allowed Topic classifier so we can upload to PyPi

* Ignore distribution egg-info build files


## 3.0.0.dev2 (4 September 2018)

### Features

* Export: Add safety checks for optional annotations and geo data

* Include more lat/longs in the default geo data

### Development

* Add release tooling

* Document the release process and a few development practices

* Travis CI: Switch to rebuilding the Docker image only for new releases

* Remove ebola, lassa, tb, WNV, and zika builds now in their own repos.  These
  builds are now available at URLs like <https://github.com/nextstrain/ebola>,
  for example.


## 3.0.0.dev1 (unreleased)

### Development

* Start versioning augur beginning with 3.0.0.  A new `augur version` command
  reports the running version.
