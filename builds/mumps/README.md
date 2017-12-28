## Mumps build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna

Should populate `fauna/data/mumps.fasta`.

#### 3. Prepare analyses
```
python mumps.prepare.py --geo global
```
This creates the file `prepared/mumps_global.json`.

```
python mumps.prepare.py --geo na
```
This creates the file `prepared/mumps_na.json`.

#### 4. Run builds
```
python mumps.process.py --json prepared/mumps_global.json --geo global
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/` for both global build.

```
python mumps.process.py --json prepared/mumps_na.json --geo na
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/` for both global build.

#### 5. Copy JSONs to auspice
```
cp auspice/mumps_* ../../../auspice/data/
```

#### 6. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
