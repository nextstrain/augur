# Adding Custom Trait Colors

Auspice uses a number of default color schemes to color the tree using meta data or values that the augur pipeline computed. In some cases these defaults are not suitable for particular type of data, and you'd like to use your own color schemes.

To specify a mapping between discrete trait values and colors, you can pass a tab-delimited file to `augur export`/`augur export v1`/`augur export v2` using `--colors`.

> _Note that it's not currently possible to specify color schemes for nucleotides, amino acids, or continuous data._

### Format

The file should contain 3 columns separated by tab characters. The first columns specifies the trait (e.g. country, host, ...), the second columns specifies the value of the trait, and the third the associated color in hexadecimal notation.

Example:
```
country	thailand	#511EA8
country	vietnam	#4928B4
country	singapore	#4334BF
country	french polynesia	#4041C7
country	american samoa	#3F50CC
country	fiji	#3F5ED0
```

The hexadecimal notation starts with a #, followed by two letters coding for red, green, blue values. If you need advice or inspiration for colors, have a look at [colorbrewer](https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3).
Pay particular not to mix tabs and spaces for column separation.
The line
```
country	american samoa	#3F50CC
```
for example has a tab between "country" and "american samoa" while "american" and "samoa" are separated by a space (because that happens to be the how this place is called).
"samoa" and "#3F50CC" are again separted by a tab.

Augur matches geographic location in lower case, so the casing doesn't matter.
Beyond that, make sure your trait (e.g. 'country') and the values (e.g. 'thailand') exactly match what's in your metadata!

### Ordering

Another reason you may want to use colors is to preserve ordering of your values. For example, if you have a trait `age_range` with values `<18`, `18-30`, `31-65` and `>65`, you'd probably like them to appear in this order in the filter and legend.

By default, auspice would sort these alphabetically, so you'll get them in a strange order. However, if you specify the order by assigning colors, the order will be preserved in auspice:

```
age_range <18   #511EA8
age_range 18-30 #4928B4
age_range 31-65 #4334BF
age_range >65   #4041C7

```

