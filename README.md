# Gelex

[![GitHub issues](https://img.shields.io/github/issues/r1cheu/gelexy?color=green)](https://github.com/r1cheu/gelexy/issues/new)

Gelex is a high-performance C++ genomic analysis toolkit designed for large-scale genomic prediction and association studies. It integrates advanced Bayesian models (BayesAlphabet series) and frequentist methods (GBLUP), optimized for memory efficiency and computational speed on massive datasets.

> [!IMPORTANT]
> This project is under active development. APIs and features are subject to change.

## 🌟 Key Features

- **Comprehensive Model Support**:
  - **Bayesian Series**: BayesA, B, C, R, RR and their counterparts for dominance effects (d) and Pi estimation (pi) (14 prior strategies in total).
  - **Frequentist Methods**: GBLUP model with REML estimation.
- **Full-stack Analysis Workflow**:
  - **Model Fitting (`fit`)**: Efficient MCMC sampling or REML execution.
  - **Genomic Prediction (`predict`)**: Generate predictions for new samples based on trained effect sizes.
  - **Association Testing (`assoc`)**: Mixed linear model GWAS with LOCO (Leave-One-Chromosome-Out) support.
  - **GRM Computation (`grm`)**: Multiple algorithms (Yang, Zeng, Vitezica) for computing Genomic Relationship Matrices.
  - **Phenotype Simulation (`simulate`)**: Simulate complex additive and dominance genetic architectures based on real genotypes.
- **Exceptional Performance**:
  - **Vectorized Acceleration**: Optimized data I/O using AVX512/AVX2 instruction sets.
  - **Multi-threaded Parallelism**: OpenMP-based parallelization across all modules.
  - **Memory Efficiency**: Memory-mapped BED file reading via `mio` with chunk-based processing support.
  - **Modern Backend**: Powered by the Eigen linear algebra library with MKL or OpenBLAS support.

## 🚀 Installation

### Prerequisites

- **CMake**: 3.20+
- **Compiler**: C++23 compatible compiler (GCC 11+, Clang 14+, or MSVC 2022+)
- **Package Management**: [pixi](https://pixi.sh) (Recommended, handles all dependencies automatically)

### Quick Install

```bash
git clone --recurse-submodules https://github.com/r1cheu/gelex.git
cd gelex
pixi run install-release
```

This will install the `gelex` binary to `~/.local/bin`.

## 💡 Quick Start

### 1. Fit a Bayesian Model (BayesR)

```bash
gelex fit \
  --bfile data/genotypes \
  --pheno data/phenotypes.tsv \
  --method R \
  --iters 10000 \
  --burnin 2000 \
  --o result/my_analysis
```

### 2. Perform Genomic Prediction

```bash
gelex predict \
  --bfile data/new_samples \
  --snp-eff result/my_analysis.snp.eff \
  --o result/predictions
```

### 3. Compute GRM (LOCO Mode)

```bash
gelex grm \
  --bfile data/genotypes \
  --method yang \
  --loco \
  --o result/grm_loco
```

### 4. Run GWAS Association Test

```bash
gelex assoc \
  --bfile data/genotypes \
  --pheno data/phenotypes.tsv \
  --grm result/grm_loco \
  --loco \
  --o result/gwas_results
```

## 🛠 Development & Testing

Run the full test suite:

```bash
pixi run test
```

For detailed development guidelines (coding style, architecture), please refer to [CLAUDE.md](CLAUDE.md).

## 📄 License

This project is licensed under the [Apache-2.0 LICENSE](LICENSE).

## 📧 Contact & Citation

If you use Gelex in your research, please cite it using the following BibTeX entry:

```bibtex
@misc{gelex,
   author = {RuLei Chen},
   year = {2026},
   note = {https://github.com/r1cheu/gelex},
   title = {Gelex: A high-performance C++ genomic analysis toolkit}
}
```

## 🙏 Acknowledgements

Gelex utilizes several libraries:

- **[Eigen](https://eigen.tuxfamily.org/)**: Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.
- **[barkeep](https://github.com/proclab/barkeep)**: Small, single C++ header to display async animations, counters, and progress bars.
- **[argparse](https://github.com/p-ranav/argparse)**: Argument parser for modern C++.
- **[mio](https://github.com/mandreyel/mio)**: An easy to use header-only cross-platform C++11 memory mapping library with an MIT license.
