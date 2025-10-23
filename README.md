# Gelex

[![GitHub issues](https://img.shields.io/github/issues/r1cheu/gelexy?color=green)](https://github.com/r1cheu/gelexy/issues/new)

Gelex is a high-performance C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides efficient computation for genomic selection with memory-mapped genotype data support.

> [!NOTE]
> This project is undergoing rapid development.

## Features

- **Bayesian Models**: BayesA, BayesB(pi), BayesC(pi), BayesRR implementations
- **Frequentist Models**: GBLUP with REML estimation
- **High Performance**: Multi-threaded computation with BLAS/LAPACK backends (MKL or OpenBLAS)
- **Memory Efficient**: Memory-mapped BED file reading via mio
- **Modern C++**: Built with C++23 standards
- **CLI Interface**: Easy-to-use command-line interface

## Installation

### Prerequisites

- CMake 3.18+
- C++23 compatible compiler (GCC 11+, Clang 14+, or MSVC 2022+)
- pixi (for dependency management)
- Git (for cloning with submodules)

### From Source

```bash
git clone --recurse-submodules https://github.com/r1cheu/gelexy.git
cd gelexy
pixi install
```

> [!NOTE]
> This project uses Git submodules for Eigen and Armadillo linear algebra libraries. All other dependencies (spdlog, mio, argparse, Catch2, etc.) are automatically managed by pixi. If you forgot to use `--recurse-submodules` during clone, run:
>
> ```bash
> git submodule update --init --recursive
> ```

### Build and Install

**Debug Build (with tests):**

```bash
pixi run build-debug
pixi run install-debug
```

**Release Build (optimized):**

```bash
pixi run build-release
pixi run install-release
```

This will install the `gelex` binary to `~/.local/bin`.

## Quick Start

### CLI Usage

Gelex provides two main subcommands:

**Fit a genomic prediction model:**

```bash
gelex fit [options]
```

**Run simulations:**

```bash
gelex simulate [options]
```

Use `gelex --help` for detailed usage information.

### Example: Bayesian Model

```bash
gelex fit \
  --bfile data/genotypes.bed \
  --pheno data/phenotypes.tsv \
  --method RR \
  --iters 10000 \
  --burnin 2000 \
  --o output
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

[License information to be added]

## Citation

[Citation information to be added]
