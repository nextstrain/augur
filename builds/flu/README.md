# Seasonal flu (H3N2, H1N1pdm, B/Vic, B/Yam) build

### Status
* all lineages working
* references only defined for HA (so far)
* subsampling is always a single pass and cannot be repeated
* tree frequencies have not been tested properly
* matchClades has not been properly tested

### summary of the necessary files

| File         | Details           |
| ------------- | ------------- |
| `flu.prepare.py`    | Prepare script.      |
| `flu.process.py`    | Process script.      |
| `flu_info.py`      | Holds (a lot of) information about sequences to drop, reference genomes etc. Used by prepare.  |
| `../../fauna/flu_<SERO>_<SEG>.fasta` | (Fauna) fasta file      |
| `colors.flu.tsv` | color maps      |


### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA file via fauna and prepare analysis
```
python flu.prepare.py --lineage h3n2 --resolution 3y
```
Running this creates the file `prepared/flu_h3n2_ha_3y.json`.

#### 3. Run build
```
python flu.process.py --json prepared/flu_h3n2_ha_3y.json
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/flu_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```

### How to run (batch)

The script `run_flu.py` will batch calls to `flu.prepare.py` and `flu.process.py` to run all combinations of lineage (`h3n2`, `h1n1pdm`, `vic`, `yam`) and resolution (`2y`, `3y`, `6y`, `12y`).
```
python run_flu.py
```
