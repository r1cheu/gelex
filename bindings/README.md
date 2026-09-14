# gelexy

Python bindings for [GELEX](https://github.com/r1cheu/gelex): genotype encoding,
genomic values from trained models, matrix reading, and posterior draws for
ArviZ. Linux x86-64, Python 3.14.

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

Compute genomic values from a trained model. The `.snplut` file written by
`gelex mcmc` holds the per-marker lookup tables; the `.bed` must carry the
training markers in the same order and allele orientation (check
`bed.marker_ids` against the `SNP`/`A1`/`A2` columns of the `.snpeff` file):

```python
import numpy as np

bed = gelexy.open_bed("data")
# Optional: reuse the training samples and their order from the .id file
# (tab-separated FID IID); any list of (FID, IID) pairs works.
with open("run.id") as f:
    bed.gather([tuple(line.split()) for line in f])
design = gelexy.make_genetic_design(bed, gelexy.load_snp_luts("run.snplut"))

sparse = gelexy.CscReader("run.draws.csc")
coefficients = sparse["genetic/A/coefficients"]  # markers × draws
projection = design.projection(gelexy.GeneticMode.A)
gebv = np.column_stack(
    [gelexy.gebv(projection, coefficients[:, d].toarray().ravel())
     for d in range(coefficients.shape[1])]
)  # samples × draws
```

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

Run `pixi run publish-python` to build, test, and publish the current version
to `https://prefix.dev/gelex` using your Pixi authentication.
