## How to run augur:
* see the `README.md` files in the respective pathogen's folder (`flu`, `ebola`, e.t.c.)

## Docs:
* [Fauna](https://github.com/nextstrain/fauna)
* [Prepare](prepare.md)
* [Process](process.md)
* [Format of (auspice) output JSONs](auspice_output.md)

## Tests

Install [pytest](https://docs.pytest.org/en/latest/).

```bash
# Install with pip.
pip install pytest

# Or install with conda.
conda install pytest
```

Checkout augur and temporarily install it as a package.

```bash
git clone --recursive https://github.com/nextstrain/augur
cd augur
pip install -e .
```

Run tests.

```bash
pytest
```

Uninstall the augur package.

```bash
pip uninstall augur
```
