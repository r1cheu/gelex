# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Gelexy is a C++ library and CLI tool for genomic prediction with Bayesian (BayesAlphabet models) and frequentist (GBLUP) approaches. It provides high-performance computation for genomic selection with memory-mapped genotype data support.

**CLI Subcommands:**
- `gelex fit`: Fit genomic prediction models
- `gelex simulate`: Run simulations

## Key Architecture

- **Core C++ Library**: Located in `src/` with modern C++23 implementation
- **Data Handling**: BED file readers, GRM computation, genotype/phenotype alignment
- **Models**: Bayesian (BayesA, BayesB, BayesC, BayesRR) and GBLUP
- **Estimation**: MCMC for Bayesian models, REML for GBLUP
- **Linear Algebra**: Uses Eigen (primary) and Armadillo (legacy) with BLAS/LAPACK backends

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

Example: `pixi run configure --build_dir=.build_release --build_type=Release`

**CMake Options:**

- `-DUSE_MKL=ON/OFF`: Use Intel MKL or OpenBLAS (default: ON in pixi tasks)
- `-DBUILD_TEST=ON/OFF`: Enable testing (default: ON in pixi tasks)

**Single Test Execution:**

```bash
cd .build_debug/tests/
./test "TestName*"
```

## Key Directories

- `src/`: Core C++ implementation
- `include/gelex/`: C++ headers
- `tests/`: Test suites using Catch2
- `ext/`: External dependencies (Eigen, Armadillo, mio)

## Core Components

- **Data Processing**: `src/data/` - BED file readers, GRM computation, data pipelines
- **Bayesian Models**: `src/model/bayes/` - BayesAlphabet implementations
- **Estimation**: `src/estimator/` - MCMC and REML estimators
- **Optimization**: `src/optim/` - Optimization policies and algorithms
- **Logging**: `src/logger/` - Performance and diagnostic logging

## Data Handling

- **BedPipe**: Memory-mapped BED file reader (`src/data/bedpipe.cpp`)
- **GRM**: Genomic Relationship Matrix computation (`src/data/grm.cpp`)
- **DataPipe**: Data processing pipeline (`src/data/datapipe.cpp`)

## Model Architecture

- **BayesModel**: Main Bayesian model class (`include/gelex/model/bayes/model.h`)
- **Genetic Traits**: BayesA/B/C/RR implementations (`src/model/bayes/traits/`)
- **Effects Management**: Fixed, random, additive, dominant effects

## Performance Features

- Memory-mapped genotype data via `mio`
- Multi-threaded computations with OpenMP
- BLAS/LAPACK acceleration (MKL or OpenBLAS)
- Chunk-based processing for large datasets

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

## Dependencies

- **Eigen**: Primary linear algebra library
- **Armadillo**: Legacy linear algebra (being phased out)
- **spdlog**: Logging library
- **mio**: Memory-mapped I/O
- **argparse**: Command-line argument parsing
- **Catch2**: Testing framework

## Code Style & Standards

- **C++ Standard**: C++23
- **CMake**: 3.18+
- **Naming**: CamelCase for classes, snake_case for functions/variables
- **Memory**: Smart pointers and RAII patterns
- **Concurrency**: OpenMP for parallel computation

## Common Development Tasks

- Adding new Bayesian models: Implement in `src/model/bayes/traits/`
- Modifying MCMC: See `src/estimator/bayes/mcmc.cpp`
- Adding data readers: Follow patterns in `src/data/`
- Adding tests: Use Catch2 framework in `tests/`
