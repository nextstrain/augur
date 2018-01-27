## Dengue build

Can run either a serotype-specific build or an all-serotype build.

**N.B.:** By default, serotype-specific builds use the respective reference genomes specified by [LANL](https://hfv.lanl.gov/content/sequence/HFV/GenomeMapper/GenomeMapper.html).
The all-serotype build uses the reference genome from serotype 4.

### How to run

#### 1. Run all commands from this directory

#### 2. Download sequence FASTAs file via fauna

Requires 5 FASTAs (`dengue_all.fasta`, `dengue_denv1.fasta`, `dengue_denv2.fasta`, `dengue_denv3.fasta` and `dengue_denv4.fasta`).

#### 3. Use custom titer tsv file from Sidney:

Requires `dengue_titers.tsv`.

#### 4. Prepare analysis
```
python dengue.prepare.py
```
Running `dengue.prepare.py` creates the files `prepared/dengue_all.json`, `prepared/dengue_denv1.json`, etc...

#### 5. Run build
```
python dengue.process.py
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 6. Copy JSONs to auspice
```
cp auspice/dengue_* ../../../auspice/data/
```

#### 7. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
