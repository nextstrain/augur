## H7N9 build

### How to run

#### 1. Run all commands from this directory

#### 2. Download FASTA files for all segments via fauna and prepare analysis
```
python avian.prepare.py
```
Running `avian.prepare.py` creates the files `prepared/avian_ha.json`, `prepared/avian_na.json`, etc...

#### 3. Run builds
```
python avian.process.py --json prepared/avian_h7n9_ha.json
python avian.process.py --json prepared/avian_h7n9_na.json
python avian.process.py --json prepared/avian_h7n9_ns.json
python avian.process.py --json prepared/avian_h7n9_pb1.json
python avian.process.py --json prepared/avian_h7n9_pb2.json
python avian.process.py --json prepared/avian_h7n9_mp.json
python avian.process.py --json prepared/avian_h7n9_np.json
python avian.process.py --json prepared/avian_h7n9_pa.json
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/avian_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
