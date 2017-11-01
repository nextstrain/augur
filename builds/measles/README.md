## Mumps build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna

Should populate `fauna/data/measles.fasta`.

#### 3. Prepare analysis
```
python measles.prepare.py
```
This creates the file `prepared/measles.json`.

#### 4. Run build
```
python measles.process.py
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/` for both global build.

#### 5. Copy JSONs to auspice
```
cp auspice/measles_* ../../../auspice/data/
```

#### 6. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
