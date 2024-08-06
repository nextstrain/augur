# Adding Custom Lat-Long Data

To place sequences on a map, auspice needs to know where the locations you specify are on the globe. These coordinates are provided by augur and augur has preset coordinates for many places. For example, for `region`, `country`, and some `division`s augur already knows many lat-long coordinates (see which ones it already knows by checking the list [here](https://github.com/nextstrain/augur/blob/-/augur/data/lat_longs.tsv)).

Some places, however, are inevitably missing and you need to tell augur the latitudes and longitude of all locations not covered by defaults. The file for latitude and longitude needs to be a tab-delimited table with four columns:

```
country	brazil	-10.3333332	-53.1999999
country	colombia	2.893108	-73.7845142
country	dominican republic	18.50012	-69.98857
country	ecuador	-1.3397667	-79.3666964
country	french polynesia	-17.6797	-149.4068
country	guatemala	15.6356088	-89.8988086
[...]
```

The first columns specifies the type of location (country, city, region, ...) corresponding to a field in your meta data. The second is the actual name of the location, while the third and fourth columns are latitude and longitude.
Note that multipart country names like "Dominican Republic" can contain a space because the delimiter of difference columns are tabs.

Augur matches geographic location in lower case, so the casing doesn't matter.
Beyond that, make sure your trait (e.g. 'country') and the values (e.g. 'brazil') exactly match what's in your metadata!

You can add geographic lat-long data at as many resolutions as you like - just ensure you give each different names ('country', 'region', 'state', 'location', 'city') and that these match columns in your metadata file.

