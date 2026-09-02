# GELEX: GEnome LEX

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Documentation Status](https://readthedocs.org/projects/gelex/badge/?version=latest)](https://gelex.readthedocs.io/en/latest/?badge=latest)
[![CI](https://github.com/r1cheu/gelex/actions/workflows/ci.yml/badge.svg)](https://github.com/r1cheu/gelex/actions/workflows/ci.yml)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/r1cheu/gelex)](https://github.com/r1cheu/gelex/releases)
[![GitHub issues](https://img.shields.io/github/issues/r1cheu/gelex)](https://github.com/r1cheu/gelex/issues)

<p align="center">
  <img src="docs/images/gelex_logo.jpeg" alt="Gelex Logo" width="640">
</p>

Gelex is a C++ toolkit for genomic prediction and genome-wide association
studies. It provides Bayesian whole-genome regression (the BayesAlphabet
family), variance component estimation, mixed-model association testing, and
supporting utilities for genomic relationship matrices, all operating directly
on PLINK binary genotypes.

> [!IMPORTANT]
> This project is under active development; interfaces and outputs may change
> between releases.

## Commands

| Command    | Description                                           |
| ---------- | ----------------------------------------------------- |
| `mcmc`     | Fit genomic prediction models with MCMC               |
| `reml`     | Estimate variance components with REML (AI algorithm) |
| `assoc`    | Mixed-model association testing, with optional LOCO   |
| `predict`  | Predict phenotypes from fitted SNP effects            |
| `grm`      | Build genomic relationship matrices from PLINK data   |
| `post`     | Summarize posterior diagnostics from MCMC samples     |

Available priors, effect modes, and options for each command are described in
the [documentation](https://gelex.readthedocs.io).

## Installation

Packages are published to the `gelex` channel on [prefix.dev](https://prefix.dev);
conda-forge supplies the runtime dependencies. Using [pixi](https://pixi.sh):

```bash
pixi global install -c conda-forge -c https://prefix.dev/gelex gelex
```

Or with [conda](https://docs.conda.io/en/latest/):

```bash
conda install -c conda-forge -c https://prefix.dev/gelex gelex
```

See the [installation guide](https://gelex.readthedocs.io/en/latest/installation.html)
for building from source.

## Usage

Fitting a BayesR model:

```bash
gelex mcmc \
  --bfile data/genotypes \
  --pheno data/phenotypes.tsv \
  --method R \
  --iters 10000 \
  --burn-in 2000 \
  --out result/my_analysis
```

Tutorials covering prediction, relationship matrices, and association testing
are available in the [documentation](https://gelex.readthedocs.io).

## Development

```bash
pixi run test
```

Coding conventions and architecture notes are documented in [CLAUDE.md](CLAUDE.md).

## License

Distributed under the [Apache-2.0 License](LICENSE).

## Citation

```bibtex
@misc{gelex,
   author = {RuLei Chen},
   year = {2026},
   note = {https://github.com/r1cheu/gelex},
   title = {Gelex: A C++ toolkit for genomic prediction and association studies}
}
```

## Acknowledgements

Gelex builds on [Eigen](https://eigen.tuxfamily.org/),
[CLI11](https://github.com/CLIUtils/CLI11),
[mio](https://github.com/mandreyel/mio), and
[barkeep](https://github.com/proclab/barkeep).
