## Zika build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna and prepare analysis
```
python zika.prepare.py
```
Alternatively, supply a FASTA file with conforming headings and pass in this file with `--sequences`.
```
python zika.prepare.py --sequences example_data/zika.fasta
```
Running `zika.prepare.py` creates the file `prepared/zika.json`.

#### 3. Run build
```
python zika.process.py
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/zika_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
