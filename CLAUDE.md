# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Gelex is a C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides high-performance computation for genomic selection with memory-mapped genotype data support.

**CLI Subcommands:**
- `gelex fit`: Fit genomic prediction models (Bayesian or GBLUP)
- `gelex simulate`: Run simulations

**Note:** The `gelex predict` subcommand mentioned in README is not yet implemented in the current codebase, but prediction utilities are available in `src/predictor/`.

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

**CMake Configuration:**

The `pixi run configure` task accepts arguments:
- `build_dir`: Build directory (default: `.build_debug`)
- `build_type`: Build type (default: `Debug`, options: `Debug`, `Release`, `RelWithDebInfo`, `MinSizeRel`)
- `test`: Enable testing (default: `ON`)
- `cxx_compiler`: C++ compiler (default: `clang++-21`)

Example: `pixi run configure --build_dir=.build_release --build_type=Release`

**Note:** The pixi configuration automatically sets `USE_MKL=ON` for optimal performance.

**CMake Options:**

- `-DUSE_MKL=ON/OFF`: Use Intel MKL or OpenBLAS (default: ON in pixi tasks)
- `-DBUILD_TEST=ON/OFF`: Enable testing (default: ON in pixi tasks)

**Compiler Configuration:**
- Uses Clang 21+ as default compiler
- OpenMP support enabled for parallel computation

**Single Test Execution:**

```bash
cd .build_debug/tests/
./test "TestName*"
```

**Test Listing:**

```bash
cd .build_debug/tests/
./test --list-tests
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
- **Prediction**: `src/predictor/` - SNP effect processing and covariate handling

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

- **DataPipe**: Main data processing pipeline (`include/gelex/data/data_pipe.h`) - handles phenotype, covariates, and sample alignment
- **BedPipe**: Memory-mapped BED file reader (`src/data/bedpipe.cpp`)
- **GenotypePipe**: Genotype data processing pipeline (`src/data/genotype_pipe.cpp`)
- **GRM**: Genomic Relationship Matrix computation (`src/data/grm.cpp`)
- **SampleManager**: Sample ID management and alignment (`src/data/sample_manager.cpp`)

## Model Architecture

- **BayesModel**: Main Bayesian model class (`include/gelex/model/bayes/model.h`) - manages fixed, random, additive, and dominant effects
- **MCMC**: Markov Chain Monte Carlo sampler (`include/gelex/estimator/bayes/mcmc.h`) - template-based implementation for different trait samplers with multi-chain parallel execution
- **Trait Samplers**: Additive effect samplers (`src/model/bayes/samplers/additive/`) - A, B, C, R, RR, RRD implementations using Gibbs and Metropolis-Hastings sampling
- **Effects Management**: Fixed, random, additive, dominant effects via `BayesEffects` system (`src/types/bayes_effects.h`)
- **BayesState**: State management for MCMC sampling (`include/gelex/model/bayes/model.h`)

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
3. **Test**: `pixi run test` or run individual tests from `.build_debug/tests/`
4. **Build Release**: `pixi run build-release` (for production builds)
5. **Install**: `pixi run install-debug` or `pixi run install-release`

## Testing

- **C++ Tests**: Catch2 framework in `tests/` directory
- **Test Data**: Uses sample BED files for validation
- **Test Execution**:
  - All tests: `pixi run test`
  - Single test: `cd .build_debug/tests && ./test "TestName*"`
  - List tests: `cd .build_debug/tests && ./test --list-tests`

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
- **Naming**: CamelCase for classes, snake_case for functions/variables
- **Memory**: Smart pointers and RAII patterns
- **Concurrency**: OpenMP for parallel computation

## Common Development Tasks

- **Adding new Bayesian models**: Implement trait samplers in `src/model/bayes/samplers/additive/`
- **Modifying MCMC**: See `include/gelex/estimator/bayes/mcmc.h` and `src/estimator/bayes/mcmc.cpp`
- **Adding data readers**: Follow patterns in `src/data/` (use `DataPipe` for data processing)
- **Adding CLI options**: Modify `src/cli/fit_command.cpp` for new command-line arguments
- **Adding tests**: Use Catch2 framework in `tests/` with sample BED files
- **Performance optimization**: Use chunk-based processing and memory-mapped I/O for large datasets
- **Adding prediction features**: Extend `src/predictor/` components for new prediction capabilities

## Important Implementation Details

- **Memory Management**: The project uses RAII patterns with smart pointers for automatic resource management
- **Error Handling**: Uses `std::expected` for error propagation with custom `Error` type
- **Parallelism**: Multi-threaded MCMC with OpenMP, with explicit thread management in MCMC implementation
- **Template-based Design**: MCMC is template-based to support different trait samplers
- **Header Organization**: Public headers in `include/gelex/`, implementation details in `src/`
- **Effect System**: Hierarchical effect management with `FixedEffect`, `RandomEffect`, `AdditiveEffect`, `DominantEffect`, and `Residual` components
- **State Tracking**: `BayesState` maintains current parameter states during MCMC sampling
- **Convergence Diagnostics**: Built-in convergence monitoring and diagnostic tools
- **Data Alignment**: `SampleManager` ensures consistent sample ordering across genotype and phenotype data
- **Chunk Processing**: Large genotype datasets processed in configurable chunks to control memory usage
