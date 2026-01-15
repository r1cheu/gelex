# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Test Commands

```bash
# Setup and build
pixi install                     # Install dependencies
pixi run build-debug            # Debug build with tests (default)
pixi run build-release          # Optimized release build
pixi run build-native           # Release with -march=native

# Testing
pixi run test                   # Run all C++ tests (ctest)
./build/debug/tests/gelex_tests --list-tests   # List all tests
./gelex_tests "TestName*"       # Run tests matching pattern

# Installation
pixi run install-debug          # Install debug binary to ~/.local/bin
pixi run install-release
```

## Code Conventions

- **C++ Standard**: C++23
- **Naming**: `PascalCase` for classes, `snake_case` for functions/variables
- **Private members**: Trailing underscore (`variable_`)
- **Function syntax**: Use trailing return type (`auto func() -> void`)
- **Include guards**: `#ifndef ... #define ... #endif`
- **Import order**: 1) Standard library, 2) External (Eigen, spdlog), 3) Internal headers
- **Eigen references**: Use `Eigen::Ref<EigenType>` and `const Eigen::Ref<const EigenType>&`
- **Comments**: Only when code is non-standard or hard to understand
- **Error handling**: Use exceptions from `include/gelex/exception.h`
- **Memory**: RAII patterns with smart pointers; use memory-mapped I/O (mio) for large files

## Architecture Overview

### Core Modules

- **CLI** (`src/cli/`): Entry points (`fit_command`, `simulate_command`, `grm_command`)
- **Data** (`src/data/`): BED file reading (`BedPipe`), genotype loading, sample alignment, GRM computation
- **Model** (`src/model/`): `BayesModel` (BayesA/B/C/R/RR/RRD), `GBLUP` (frequentist)
- **Estimator** (`src/estimator/`): MCMC samplers for Bayesian models, EM optimizer for GBLUP
- **Predict** (`src/predict/`): Prediction engine with genotype alignment and SNP matching

### Key Design Patterns

- **Pipeline**: `DataPipe`, `BedPipe` for chunk-based data processing
- **Strategy**: `VariantProcessor`, `CodePolicy` for GRM computation
- **Factory**: `create_prior_strategy()` for BayesAlphabet prior selection

### Data Flow

```
BED/BIM/FAM → BedPipe → GenotypeLoader → Model → Estimator → Output
Phenotype/Covariate → SampleManager (alignment) → Model
```

## Key Dependencies

- **Linear algebra**: Eigen with MKL/OpenBLAS backend
- **Testing**: Catch2 v3
- **Memory mapping**: mio
- **CLI**: argparse
- **Parallelism**: OpenMP

## Notes

- Public API in `include/gelex/`, implementation in `src/`
- PLINK BED format for genotypes, TSV for phenotypes
- BayesAlphabet models: A, B, Bpi, C, Cpi, R, RR, RRD (plus dominance variants)
