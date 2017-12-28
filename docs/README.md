## How to run augur:
* see the `README.md` files in the respective pathogen's folder (`flu`, `ebola`, e.t.c.)

## Directory structure

Throughout we assume that the fauna, augur and auspice are sister directories.

## Docs:
* [Fauna](https://github.com/nextstrain/fauna)
* [Prepare](prepare.md)
* [Process](process.md)
* [Format of (auspice) output JSONs](auspice_output.md)

## Tests

Tests run using [tox](http://tox.readthedocs.io/en/latest/)
and [pytest](https://docs.pytest.org/en/latest/). If you didn't install augur
with `pip -r requirements.txt` or with
the
[janus conda environment](https://github.com/nextstrain/janus/#installation),
install them as follows.

```bash
# Install with pip.
pip install --user pytest tox

# Or install inside a conda environment.
source activate janus_python2
conda install pytest virtualenv
pip install tox
```

Run tests from the augur root directory.

```bash
tox
```

Tox builds an empty virtual environment, installs augur's requirements and a
package of augur built from the local directory, and runs tests with pytest. The
`tox.ini` file defines the test environment and the `pytest.ini` file defines
which tests to run. Tox reuses the virtual environment with subsequent runs. Use
the `-r` flag to force it to rebuild its environment.

```bash
tox -r
```

Note that on OS X, matplotlib may throw errors during testing stating

```bash
from matplotlib.backends import _macosx
E   RuntimeError: Python is not installed as a framework.
```

This can be resolved by switching matplotlib backends. Create a file in the home directory `~/.matplotlib/matplotlibrc` and include the following:

```python
backend: TkAgg
```

This should fix the error.
