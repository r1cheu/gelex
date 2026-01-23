# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Install dependencies
pixi install

# Build
pixi run build-debug           # Debug build with tests and sanitizers
pixi run build-release         # Release build
pixi run build-native          # Release with -march=native optimization

# Test
pixi run test                  # Run all tests
./build/debug/tests/gelex_tests "TestName*"   # Run specific test

# Install
pixi run install-release       # Install to ~/.local/bin
```

## Architecture

**Gelex** is a high-performance genomic prediction library implementing Bayesian (BayesA/B/C/R/RR/RRD) and Frequentist (GBLUP) methods.

### Core Layers

- **Data Layer** (`src/data/`): Loaders for PLINK files (BED/BIM/FAM), phenotypes, covariates; GRM computation with LOCO support; memory-mapped genotype access
- **Model Layer** (`src/model/`): Bayesian models with Gibbs/MH samplers (`src/model/bayes/samplers/`); GBLUP model (`src/model/freq/`)
- **Optimizer Layer** (`src/optim/`): REML convergence, variance component estimation
- **Estimator Layer** (`src/estimator/`): MCMC runners, posterior calculation, REML solver
- **GWAS Layer** (`src/gwas/`): Association testing with optional LOCO correction
- **Predict Layer** (`src/predict/`): Genomic prediction engine, SNP matching and allele alignment

### CLI Commands

The CLI (`apps/`) provides five subcommands:
- `fit` - Fit Bayesian/Frequentist models
- `predict` - Genomic prediction using trained models
- `grm` - Compute genomic relationship matrix
- `assoc` - GWAS association testing
- `simulation` - Simulate genomic data

Each command has separate `*_command.{h,cpp}` (execution) and `*_args.{h,cpp}` (argument parsing) files.

### Key Dependencies

- **Eigen** (`ext/eigen/`): Linear algebra (submodule)
- **OpenMP**: Parallelization
- **BLAS/LAPACK**: MKL or OpenBLAS backend
- **spdlog**: Logging
- **argparse**: CLI parsing (`ext/argparse.h`)
- **Catch2**: Testing framework

## Code Style

- **Classes**: `PascalCase`
- **Functions/variables**: `snake_case`
- **Private members**: trailing underscore (e.g., `variable_`)
- **C++ standard**: C++23
- **Return types**: Use trailing syntax `auto func() -> ReturnType`
- **Comments**: Only when necessary (non-obvious logic)
- **Format**: Chromium style, 80 columns, Allman braces
