# Augur manuscript

## Setup

``` bash
conda create -n augur-manuscript graphviz pandoc
conda activate augur-manuscript
```

## Compile manuscript

Use pandoc to compile the paper.

```bash
make
```

Alternately, use [JOSS's Docker image to compile the paper](https://joss.readthedocs.io/en/latest/submitting.html#docker) from this `paper` directory.

``` bash
docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft
```
