# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Reference

**Common Commands:**

```bash
# Build and test
pixi run build-debug    # Debug build with tests
pixi run test           # Run all tests
pixi run build-release  # Optimized release build

# Installation
pixi run install-debug   # Install debug binary
pixi run install-release # Install release binary

# Testing
cd build/debug/tests/
./gelex_tests "TestName*"       # Run specific test
./gelex_tests --list-tests      # List all tests
```

## Project Overview

Gelex is a C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides high-performance computation for genomic selection with memory-mapped genotype data support.

**CLI Subcommands:**

- `gelex fit`: Fit genomic prediction models (Bayesian or GBLUP)
- `gelex simulate`: Run simulations

**Note:** The `gelex predict` subcommand mentioned in README is not yet implemented in the current codebase, but prediction utilities are available in `src/predict/`.

## Key Architecture

- **Core C++ Library**: Located in `src/` with modern C++23 implementation
- **Data Handling**: BED file readers, GRM computation, genotype/phenotype alignment via `DataPipe`
- **Models**: Bayesian (BayesA, BayesB, BayesC, BayesR, BayesRR, BayesRRD) and GBLUP
- **Estimation**: MCMC for Bayesian models, REML for GBLUP
- **Linear Algebra**: Uses Eigen (primary) and Armadillo (legacy) with BLAS/LAPACK backends
- **Memory Management**: Memory-mapped genotype data with chunk-based processing
- **Parallelism**: Multi-threaded MCMC with OpenMP, vectorized genotype processing

## Build System

This project uses CMake with pixi for dependency management:

**Common Commands:**

- Debug build: `pixi run build-debug` (default, includes tests)
- Release build: `pixi run build-release` (optimized)
- Install debug: `pixi run install-debug` (copies to ~/.local/bin)
- Install release: `pixi run install-release` (copies to ~/.local/bin)
- Test: `pixi run test` (runs all C++ tests with ctest)

**CMake Presets:**

The project uses CMake presets for different build configurations:

- `debug`: Debug build with tests and sanitizers enabled (default)
- `release`: Optimized release build
- `release-native`: Release build with `-march=native` optimization
- `coverage`: Debug build with code coverage instrumentation

**Note:** The pixi configuration automatically sets `USE_MKL=ON` for optimal performance.

**Compiler Configuration:**

- Uses Clang 21+ as default compiler
- OpenMP support enabled for parallel computation

**Test Execution:**

```bash
# Run all tests
pixi run test

# Run single test
cd build/debug/tests/
./gelex_tests "TestName*"

# List all tests
cd build/debug/tests/
./gelex_tests --list-tests

# Generate coverage report
pixi run coverage
```

## Key Directories

- `src/`: Core C++ implementation
- `include/gelex/`: C++ headers
- `tests/`: Test suites using Catch2
- `ext/`: External dependencies (Eigen, Armadillo, mio)

## Core Components

- **Data Processing**: `src/data/` - BED file readers, GRM computation, data pipelines (`DataPipe`, `BedPipe`, `GenotypePipe`, `SampleManager`)
- **Bayesian Models**: `src/model/bayes/` - BayesAlphabet implementations (`BayesModel`, trait samplers, `BayesState`)
- **Frequentist Models**: `src/model/freq/` - GBLUP implementation
- **Estimation**: `src/estimator/` - MCMC for Bayesian models, REML for GBLUP
- **Optimization**: `src/optim/` - Optimization policies and algorithms
- **Logging**: `src/logger/` - Performance and diagnostic logging
- **CLI Interface**: `src/cli/` - Command-line interface implementation
- **Utilities**: `src/utils/` - Math utilities, formatters, and helper functions
- **Prediction**: `src/predict/` - SNP effect processing and covariate handling

## Data Flow Architecture

1. **Data Loading**: `DataPipe` loads phenotype/covariates, `BedPipe` loads genotypes via memory-mapping
2. **Sample Alignment**: `SampleManager` aligns samples across all data sources
3. **Model Setup**: `BayesModel` initializes effects (fixed, random, additive, dominant)
4. **Estimation**: `MCMC` samples from posterior distribution or `REML` estimates variance components
5. **Output**: Results written to disk with convergence diagnostics

**Key Data Structures**:

- `BayesEffects`: Manages fixed, random, additive, dominant, and residual effects
- `BayesState`: Tracks current state of all model parameters during MCMC
- `MCMCResult`: Stores posterior samples and convergence statistics

## Data Handling

- **DataPipe**: Main data processing pipeline - handles phenotype, covariates, and sample alignment
- **BedPipe**: Memory-mapped BED file reader
- **GenotypePipe**: Genotype data processing pipeline
- **GRM**: Genomic Relationship Matrix computation
- **SampleManager**: Sample ID management and alignment

## Model Architecture

- **BayesModel**: Main Bayesian model class - manages fixed, random, additive, and dominant effects
- **MCMC**: Markov Chain Monte Carlo sampler - template-based implementation for different trait samplers with multi-chain parallel execution
- **Trait Samplers**: Additive effect samplers - A, B, C, R, RR, RRD implementations using Gibbs and Metropolis-Hastings sampling
- **Effects Management**: Fixed, random, additive, dominant effects via `BayesEffects` system
- **BayesState**: State management for MCMC sampling

**Bayesian Model Types**:

- **BayesA**: t-distributed marker effects
- **BayesB**: Mixture model with spike-slab prior
- **BayesBpi**: Bayes B with variable π
- **BayesC**: Multi-component mixture
- **BayesCpi**: Bayes C with variable π
- **BayesR**: Multi-component with fixed variances
- **BayesRR**: Ridge regression (Bayesian LASSO)
- **BayesRRD**: Ridge regression with dominance effects

## Performance Features

- **Memory-mapped genotype data** via `mio` library with AVX512/AVX2 vectorization
- **Multi-threaded computations** with OpenMP for MCMC chains and matrix operations
- **BLAS/LAPACK acceleration** (MKL or OpenBLAS) for linear algebra
- **Chunk-based processing** for large datasets to control memory usage
- **Vectorized genotype processing** with automatic monomorphic variant detection
- **Cross-GRM computation** for prediction scenarios

## Development Workflow

1. **Setup**: `pixi install` (installs all dependencies)
2. **Build Debug**: `pixi run build-debug` (includes configure step)
3. **Test**: `pixi run test` or run individual tests from `build/debug/tests/`
4. **Build Release**: `pixi run build-release` (for production builds)
5. **Install**: `pixi run install-debug` or `pixi run install-release`

**Debug Build Features:**

- AddressSanitizer and UndefinedBehaviorSanitizer enabled
- Full debug symbols
- Test binaries included in build directory

## Testing

- **C++ Tests**: Catch2 framework in `tests/` directory
- **Test Data**: Uses sample BED files for validation
- **Test Execution**:
  - All tests: `pixi run test`
  - Single test: `cd build/debug/tests && ./gelex_tests "TestName*"`
  - List tests: `cd build/debug/tests && ./gelex_tests --list-tests`
  - Run by tag: `cd build/debug/tests && ./gelex_tests [tag]`

**Test Organization**:

- `test_*.cpp` files cover specific components (data loading, GRM computation, MCMC, etc.)
- Tests use sample PLINK binary files for realistic validation
- Each test module focuses on a specific subsystem (data, model, estimator, etc.)

## Dependencies

- **Eigen**: Primary linear algebra library
- **Armadillo**: Legacy linear algebra (being phased out)
- **spdlog**: Logging library
- **mio**: Memory-mapped I/O
- **argparse**: Command-line argument parsing
- **Catch2**: Testing framework

## CLI Interface

- **Main Entry**: `src/main.cpp` handles command-line parsing and subcommand dispatch
- **Subcommands**: Implemented in `src/cli/` directory
  - `fit_command.cpp`: Model fitting with Bayesian and GBLUP methods
  - `simulation_command.cpp`: Simulation functionality
- **Argument Parsing**: Uses argparse library with comprehensive option validation
- **Configuration**: CLI options mapped to internal configuration objects for model setup

## Code Style & Standards

- **C++ Standard**: C++23
- **CMake**: 3.18+
- **Naming**: PascalCase for classes, snake_case for functions/variables
- **Memory**: Smart pointers and RAII patterns
- **Concurrency**: OpenMP for parallel computation

## Common Development Tasks

- **Adding new Bayesian models**: Implement trait samplers in `src/model/bayes/samplers/additive/`
- **Modifying MCMC**: See `include/gelex/estimator/bayes/mcmc.h` and `src/estimator/bayes/mcmc.cpp`
- **Adding data readers**: Follow patterns in `src/data/` (use `DataPipe` for data processing)
- **Adding CLI options**: Modify `src/cli/fit_command.cpp` for new command-line arguments
- **Adding tests**: Use Catch2 framework in `tests/` with sample BED files
- **Performance optimization**: Use chunk-based processing and memory-mapped I/O for large datasets
- **Adding prediction features**: Extend `src/predict/` components for new prediction capabilities

## Important Implementation Details

- **Memory Management**: The project uses RAII patterns with smart pointers for automatic resource management
- **Error Handling**: try-catch system in include/gelex/exception.h for robust error reporting
- **Parallelism**: Multi-threaded MCMC with OpenMP, with explicit thread management in MCMC implementation
- **Template-based Design**: MCMC is template-based to support different trait samplers
- **Header Organization**: Public headers in `include/gelex/`, implementation details in `src/`
- **Effect System**: Hierarchical effect management with `FixedEffect`, `RandomEffect`, `AdditiveEffect`, `DominantEffect`, and `Residual` components
- **State Tracking**: `BayesState` maintains current parameter states during MCMC sampling
- **Convergence Diagnostics**: Built-in convergence monitoring and diagnostic tools
- **Data Alignment**: `SampleManager` ensures consistent sample ordering across genotype and phenotype data
- **Chunk Processing**: Large genotype datasets processed in configurable chunks to control memory usage
