## Docs:
* [Prepare](prepare.md)
* [Process](process.md)
* [Format of (auspice) output JSONs](auspice_output.md)

## How to run:

### H7N9

```
# ensure H7N9 fasta files in fauna folder (see config in prepare.H7N9.py for paths)
cd H7N9
python prepare.H7N9.py # creates JSONs in the folder prepared
python process.H7N9.py # creates files in processed/ and auspice/
cp auspice/flu_H7N9_HA_* ../../auspice/data/ #assuming your directory structure is like mine!
cd ../../auspice
node dev-server.js local #or npm run start:local
```

### zika

```
# ensure zika fasta files in fauna folder (see config in zika.prepare.py for paths)
cd zika
python zika.prepare.py # creates JSONs in the folder prepared
python zika.process.py # creates files in processed/ and auspice/
cp auspice/* ../../auspice/data/ #assuming your directory structure is like mine!
cd ../../auspice
node dev-server.js local #or npm run start:local
```
