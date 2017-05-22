### instructions

* First download the data from fauna (instructions in that repo)
```
python prepare.H7N9.py
# in one terminal window:
for i in "PB1" "PB2" "PA" "NS"; do python process.H7N9.py -j prepared/flu_H7N9_${i}.json; done
# in another:
for i in "NA" "NP" "MP" "HA"; do python process.H7N9.py -j prepared/flu_H7N9_${i}.json; done
```

```
cp auspice/*json ../../auspice/data/
cd ../../auspice
node dev-server.js local
```

* to move to S3 bucket:
```
cd augur/build
rm *
cd augur/H7N9/auspice
# make everything lower case
for i in *json; do j=$( echo "$i" | tr '[:upper:]' '[:lower:]' ); cp $i ../../build/$j
bundle exec s3_website push --site build
```
