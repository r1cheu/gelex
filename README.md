# Gelex

[![GitHub issues](https://img.shields.io/github/issues/r1cheu/gelexy?color=green)](https://github.com/r1cheu/gelexy/issues/new)

Gelex is a high-performance C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides efficient computation for genomic selection with memory-mapped genotype data support.

> [!NOTE]
> This project is undergoing rapid development.

## Features

- **Bayesian Models**: BayesA, BayesB, BayesBpi, BayesC, BayesCpi, BayesR, BayesRR, BayesRRD implementations
- **Frequentist Models**: GBLUP with REML estimation
- **High Performance**: Multi-threaded computation with OpenMP, BLAS/LAPACK backends (MKL or OpenBLAS)
- **Memory Efficient**: Memory-mapped BED file reading via mio with AVX512/AVX2 vectorization
- **Modern C++**: Built with C++23 standards using RAII and smart pointers
- **Data Handling**: Chunk-based processing for large datasets, sample alignment across genotype/phenotype data
- **CLI Interface**: Easy-to-use command-line interface with comprehensive validation

## Installation

### Prerequisites

- CMake 3.20+
- C++23 compatible compiler (GCC 11+, Clang 14+, or MSVC 2022+)
- pixi (for dependency management)
- Git (for cloning with submodules)

### From Source

```bash
git clone --recurse-submodules https://github.com/r1cheu/gelex.git
cd gelex
pixi r install-release
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

### Testing

To verify your installation and run the test suite:

```bash
# Run all tests
pixi run test

# Run individual tests
cd build/debug/tests/
./gelex_tests --list-tests      # List all tests
./gelex_tests "TestName*"       # Run specific test
```

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

### Example: GBLUP Model

```bash
gelex fit \
  --bfile data/genotypes.bed \
  --pheno data/phenotypes.tsv \
  --method GBLUP \
  --o output
```

> **Note:** The `gelex predict` subcommand mentioned in some documentation is not yet implemented in the current codebase, but prediction utilities are available in the library for programmatic use.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

For development guidelines, build instructions, and codebase architecture, see [CLAUDE.md](CLAUDE.md).

## License

BSD-style license. See [LICENSE](LICENSE) file for details.

## Citation

If you use Gelex in your research, please cite:

```
Citation information to be added.
```
