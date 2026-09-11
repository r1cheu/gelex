# gelexy

Python bindings for [GELEX](https://github.com/r1cheu/gelex): genotype encoding,
matrix reading, and posterior draws for ArviZ. Linux x86-64, Python 3.14.

## Install

```bash
conda install -c https://prefix.dev/gelex -c conda-forge gelexy
```

## Usage

Read MCMC draws into ArviZ:

```python
import arviz as az
import gelexy

idata = gelexy.read_draws("run.draws")
az.summary(idata, var_names=["residual.variance"])
```

Marker-sized payloads are skipped by default. Set `include_markers=True` to
include them; sparse draws from the `.draws.csc` companion are then densified.
`DenseReader` and `CscReader` provide direct NumPy and SciPy matrix access.

Encode genotypes in place (samples × markers, writable Fortran-order `float64`):

```python
import numpy as np

x = np.array([[0.0], [1.0], [2.0]], order="F")
gelexy.encode_inplace(x, gelexy.GeneticMode.A, gelexy.GenotypeMethod.Center)
```

## Build

From the repository root:

```bash
pixi run test-python
pixi run -e package build-conda-python
```

The independent package version is in [VERSION](VERSION). Built packages are
written to `/tmp/rb-python-output/linux-64/`.
