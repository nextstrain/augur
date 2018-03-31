## Ebola build

### How to run

* Run all commands from this directory

```
python lassa.prepare.py
```
* Running `lassa.prepare.py` creates the files `prepared/lassa_S.json` and `prepared/lassa_L.json`.

#### 3. Run build
```
python lassa.process.py --json prepared/lassa_S.json.json
python lassa.process.py --json prepared/lassa_L.json.json
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
