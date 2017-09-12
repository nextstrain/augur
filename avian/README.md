### instructions

* First download the data from fauna (instructions in that repo)
* Second, run augur

```
# ensure H7N9 fasta files in fauna folder (see config in prepare.H7N9.py for paths)
cd H7N9
rm auspice/* prepared/* processed/* #makes it cleaner
python prepare.H7N9.py # creates JSONs in the folder prepared
python process.H7N9.py -j prepared/<JSON> # creates files in processed/ and auspice/
```

* As there are 8 scripts, and therefore 8 independent process scripts to run, if you have 2 cores this speeds things up:
```
# in one terminal window:
for i in "pb1" "pb2" "pa" "ns"; do python process.H7N9.py -j prepared/flu_h7n9_${i}.json; done
# in another:
for i in "na" "np" "mp" "ha"; do python process.H7N9.py -j prepared/flu_h7n9_${i}.json; done
```


* to visualise locally in auspice
```
cd augur/H7N9
cp auspice/*json ../../auspice/data/
# local build needs to be running...
cd ../../auspice
node dev-server.js local
# visit localhost:4000/flu
```

* to move to S3 bucket:
```
cd augur/build
rm *
cd augur/H7N9/auspice
bundle exec s3_website push --site build
```
