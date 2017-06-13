### instructions

* First download the data from fauna (instructions in that repo)
* Second, run augur
```
cd augur/H7N9
rm prepared/* processed/* auspice/*
python prepare.H7N9.py
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
