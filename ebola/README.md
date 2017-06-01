### how to run

* Download FASTA file via fauna and ensure the path set in `ebola.prepare.py` is correct

```
cd ebola
rm auspice/* prepared/* processed/* #makes it cleaner
python ebola.prepare.py # creates JSONs in the folder prepared
python ebola.process.py # creates files in processed/ and auspice/
cp auspice/* ../../auspice/data/ #assuming your directory structure is like mine!
cd ../../auspice
node dev-server.js local #or npm run start:local
```
