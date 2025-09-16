# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Gelexy is a Python package for genomic prediction with C++ backend for performance. It provides both Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches for genomic selection.

## Key Architecture

- **Core C++ Library**: Located in `src/` with nanobind Python bindings in `src/python/`
- **Python Interface**: `python/gelexy/` provides high-level API
- **Data Handling**: BED file readers, GRM computation, genotype/phenotype alignment
- **Models**: Bayesian (BayesA, BayesB, BayesC, BayesRR) and GBLUP
- **Estimation**: MCMC for Bayesian models, REML for GBLUP

## Build System

This project uses scikit-build-core with CMake for building Python extensions:

**Common Commands:**

- Build: `pixi run build` or `uv build`
- Configure: `pixi run configure`
- Test C++: `pixi run cpptest`
- Test Python: `pytest tests/python/`
- Single C++ test: `cd .build/tests/cpp/ && ./test --gtest_filter="TestName*"`
- Build wheel: `pixi run wheel`

**CMake Options:**

- `-DUSE_MKL=ON/OFF`: Use Intel MKL or OpenBLAS
- `-DBUILD_TEST=ON/OFF`: Enable testing
- `-DBUILD_PYTHON=ON/OFF`: Build Python bindings

## Development Workflow

1. **Setup**: `pixi install` (uses pixi.toml for dependencies)
2. **Build**: `pixi run build`
3. **Test**: `pixi run test` (runs both Python and C++ tests)
4. **Development**: Use `pixi run configure` and `pixi run build` iteratively

## Key Directories

- `src/`: Core C++ implementation
- `include/gelex/`: C++ headers
- `python/gelexy/`: Python package
- `tests/`: Test suites (C++ in `cpp/`, Python in `python/`)
- `ext/`: External dependencies (Eigen, mio) - Note: Armadillo is being phased out in favor of Eigen

## Python API Structure

- `gelexy.make_gblup()`: Create GBLUP model factory
- `gelexy.make_bayes()`: Create Bayesian model factory
- `gelexy.make_grm()`: Compute genomic relationship matrix
- `gelexy.Estimator`: Main estimation class
- `gelexy.MCMC`: Bayesian MCMC configuration

## Testing

- **C++ Tests**: Use Catch2 framework, run with `ctest`
- **Python Tests**: Use pytest, located in `tests/python/`
- **Test Data**: Located in `tests/data/` with sample BED files

## Performance Considerations

- Uses Eigen for linear algebra (Armadillo is being phased out)
- Supports both MKL and OpenBLAS backends
- Memory-mapped genotype data support via `mio`
- Multi-threaded computations with OpenMP

## Common Development Tasks

- Adding new Bayesian models: Implement in `src/model/bayes/traits/`
- Modifying MCMC: See `src/estimator/bayes/mcmc.cpp`
- Adding Python bindings: Modify `src/python/` files
- Adding tests: Follow existing patterns in `tests/`
