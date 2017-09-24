## Ebola build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna and prepare analysis
```
python ebola.prepare.py
```
Alternatively, supply a FASTA file with conforming headings and pass in this file with `--sequences`.
```
python ebola.prepare.py --sequences example_data/ebola.fasta
```
Running `ebola.prepare.py` creates the file `prepared/ebola.json`.

#### 3. Run build
```
python ebola.process.py
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/ebola_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
