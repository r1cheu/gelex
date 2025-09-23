# Gelex Project Handover Document

## Project Overview

**Gelex** is a high-performance C++ command-line application for genomic prediction, focusing on efficient processing of PLINK binary genotype files and implementing both Bayesian and GBLUP (Genomic Best Linear Unbiased Prediction) methods for genomic selection.

### Key Features

- **Bayesian Models**: Support for multiple Bayesian methods (A, B, Bpi, C, Cpi, R, RR)
- **GBLUP**: Traditional genomic relationship matrix-based prediction
- **High Performance**: Optimized C++ implementation with efficient memory management
- **PLINK Integration**: Native support for PLINK binary format (.bed, .bim, .fam)
- **Command-line Interface**: Comprehensive CLI with argparse for easy usage
- **Memory-mapped I/O**: Efficient handling of large genotype files

## Project Structure

```
gelex/
├── include/gelex/          # C++ headers
│   ├── data/              # Data handling and I/O
│   ├── model/             # Model definitions
│   ├── optim/             # Optimization interfaces
│   └── utils/             # Utilities
├── src/                   # C++ source code
│   ├── data/              # Data processing implementations
│   ├── model/             # Model implementations
│   │   ├── bayes/         # Bayesian models
│   │   └── freq/          # Frequentist models
│   ├── estimator/         # Estimation algorithms
│   │   ├── bayes/         # Bayesian estimation
│   │   └── freq/          # Frequentist estimation
│   ├── predictor/         # Prediction utilities
│   ├── optim/             # Optimization implementations
│   ├── logger/            # Logging system
│   └── utils/             # Utility implementations
├── tests/                 # Test suite
│   ├── test_*.cpp         # Unit tests for various components
│   └── CMakeLists.txt     # Test build configuration
├── ext/                   # External dependencies
│   ├── mio/               # Memory-mapped I/O library
│   ├── eigen/             # Eigen linear algebra library
│   └── armadillo/         # Armadillo numerical library
└── docs/                  # Documentation
```

## Build System & Dependencies

### Core Dependencies

- **C++**: C++23 standard
- **BLAS/LAPACK**: MKL or OpenBLAS
- **Build Tools**: CMake >=3.18, Ninja, pixi
- **External Libraries**:
  - spdlog: Logging library
  - argparse: Command-line argument parsing
  - mio: Memory-mapped I/O
  - Eigen: Linear algebra
  - Armadillo: Numerical computations

### Build Configuration

The project uses **CMake** with **pixi** for dependency management. Key build options:

- `USE_MKL=ON/OFF`: Use Intel MKL or OpenBLAS
- `BUILD_TEST=ON/OFF`: Build tests

## Core Architecture

### Main Components

1. **Data Layer** (`include/gelex/data/`, `src/data/`)
   - PLINK BED file reader
   - Genotype and phenotype data pipelines
   - GRM (Genomic Relationship Matrix) calculations
   - Data validation and processing

2. **Model Layer** (`include/gelex/model/`, `src/model/`)
   - Bayesian model implementations (A, B, Bpi, C, Cpi, R, RR)
   - GBLUP model implementation
   - Effect management system
   - Trait and distribution handling

3. **Estimation Layer** (`src/estimator/`)
   - MCMC sampling for Bayesian methods
   - Frequentist estimation algorithms
   - Diagnostics and result processing

4. **Utility Layer** (`src/utils/`, `src/logger/`)
   - Mathematical utilities
   - Logging system with multiple loggers
   - Formatter utilities

### Main Entry Point

- `src/main.cpp`: Command-line interface with argparse
- Supports "fit" subcommand for model fitting
- Handles file I/O, configuration, and pipeline execution

## Key Models & Algorithms

### Supported Methods

- **BayesA**: Single-component Bayes A
- **BayesRR**: Ridge regression Bayes
- **BayesB**: Two-component mixture
- **BayesBpi**: Bayes B with variable π
- **BayesC**: Multi-component mixture
- **BayesCpi**: Bayes C with variable π
- **BayesR**: Multi-component with fixed variances
- **GBLUP**: Genomic best linear unbiased prediction

### Data Processing

- **PLINK Integration**: Native support for .bed, .bim, .fam files
- **Memory-mapped I/O**: Efficient handling of large genotype files
- **Chunk Processing**: Configurable chunk size for SNP processing
- **Sample Alignment**: Automatic alignment of genotype and phenotype data

## Usage Examples

### Command-line Usage

```bash
# Fit a Bayesian RR model
./gelex fit --bfile genotype_data --pheno phenotype_data.tsv --pheno-col 3 --method RR --out results

# Fit GBLUP with dominance effects
./gelex fit --bfile genotype_data --pheno phenotype_data.tsv --dom --method GBLUP --out results

# Fit with covariates
./gelex fit --bfile genotype_data --pheno phenotype_data.tsv --qcovar covariates.tsv --method Bpi --out results
```

### Supported Options

- `--bfile`: PLINK binary file prefix
- `--pheno`: Phenotype file (TSV format)
- `--pheno-col`: Phenotype column index (default: 3)
- `--method`: Prediction method (A, B, Bpi, C, Cpi, R, RR, GBLUP)
- `--dom`: Enable dominance effects
- `--qcovar`: Quantitative covariates file
- `--covar`: Categorical covariates file
- `--chunk-size`: SNP processing chunk size (default: 10000)
- `--out`: Output prefix

## Development Workflow

### Environment Setup

```bash
# Using pixi (recommended)
pixi install
pixi run build

# Manual build with CMake
cmake -B build -DUSE_MKL=ON -DBUILD_TEST=ON
cmake --build build -j16
```

### Testing

```bash
# Run tests
pixi run test

# Or manually
cd build && ctest --output-on-failure
```

### Installation

```bash
# Install to local bin directory
pixi run install
```

## Performance Considerations

### BLAS Configuration

- **MKL**: Better performance on Intel systems
- **OpenBLAS**: Good cross-platform performance
- Configure with `USE_MKL=ON/OFF` in CMake

### Memory Management

- Uses memory-mapped I/O for large genotype files
- Configurable chunk processing for memory efficiency
- Efficient matrix operations with Eigen/Armadillo

### Parallelization

- OpenMP support for multi-core operations
- Thread-safe implementations where appropriate

## Key Files & Components

### Core Implementation Files

- `src/main.cpp`: Main CLI entry point
- `src/data/`: Data processing and I/O
- `src/model/bayes/`: Bayesian model implementations
- `src/model/freq/`: Frequentist model implementations
- `src/estimator/bayes/`: Bayesian estimation algorithms
- `src/estimator/freq/`: Frequentist estimation algorithms

### Configuration Files

- `CMakeLists.txt`: C++ build configuration
- `pixi.toml`: Dependency management and task runner
- `.clang-format`: Code formatting rules
- `.clang-tidy`: Static analysis configuration

## Troubleshooting & Common Issues

### Build Issues

1. **Missing BLAS/LAPACK**: Ensure MKL or OpenBLAS is installed
2. **CMake version**: Requires CMake 3.18+
3. **Compiler support**: Requires C++23 compatible compiler

### Runtime Issues

1. **File permissions**: Ensure read access to data files
2. **Memory**: Large datasets may require significant RAM
3. **File formats**: Verify PLINK binary file compatibility

## Future Development Directions

### Immediate Priorities

1. **Documentation**: Complete user documentation and examples
2. **Testing**: Expand test coverage for all components
3. **Performance**: Optimize critical computation paths

### Medium-term Goals

1. **Additional Formats**: Support for more genotype formats
2. **GPU Acceleration**: CUDA/OpenCL support for matrix operations
3. **Pipeline Integration**: Nextflow/Snakemake support

### Long-term Vision

1. **Cloud Deployment**: Containerized deployment options
2. **Web Interface**: REST API and web dashboard
3. **Extended Models**: Support for more complex genetic models

## Contact & Support

- **Maintainer**: r1cheu (chenrulei@cemps.ac.cn)
- **Repository**: https://github.com/r1cheu/gelex
- **Issues**: GitHub issue tracker

---

_This handover document was updated on 2025-09-22. Please keep it updated as the project evolves._
