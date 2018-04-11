## Lassa build

### How to run

* Run all commands from this directory

```
python lassa.prepare.py
```
* Running `lassa.prepare.py` creates the files `prepared/lassa_s.json` and `prepared/lassa_l.json`.

#### 3. Run build
```
python lassa.process.py --json prepared/lassa_s.json
python lassa.process.py --json prepared/lassa_l.json
```
This creates intermediary files in `processed/` and auspice-ready JSONs in `auspice/`.

#### 4. Copy JSONs to auspice
```
cp auspice/lassa_* ../../../auspice/data/
```

#### 5. Run auspice to visualize
```
cd ../../../auspice
npm run start:local
```
