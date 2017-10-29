## Mumps build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna and prepare analysis
```
python mumps.prepare.py
```
Running `mumps.prepare.py` creates the file `prepared/mumps.json`.

#### 3. Run build
```
python mumps.process.py
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/mumps_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
