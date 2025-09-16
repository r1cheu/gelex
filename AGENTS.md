# Gelexy Project Handover Document

## Project Overview

**Gelexy** is a high-performance Python package for genomic prediction, leveraging C++ for computational efficiency while maintaining Python usability. The project focuses on Bayesian methods and GBLUP (Genomic Best Linear Unbiased Prediction) models for genomic selection.

### Key Features

- **Bayesian Models**: Support for multiple Bayesian alphabets (A, RR, B, Bpi, C, Cpi, R)
- **GBLUP**: Traditional genomic relationship matrix-based prediction
- **High Performance**: C++ backend with Python bindings using nanobind
- **Modern Build System**: scikit-build-core with CMake integration
- **Cross-Platform**: Supports Linux with MKL or OpenBLAS backends

## Project Structure

```
gelexy/
├── include/gelex/          # C++ headers
│   ├── model/             # Model definitions
│   ├── data/              # Data handling
│   ├── optim/             # Optimization
│   └── utils/             # Utilities
├── src/                   # C++ source code
│   ├── model/             # Model implementations
│   ├── data/              # Data processing
│   ├── python/            # Python bindings
│   ├── optim/             # Optimization algorithms
│   └── logger/            # Logging utilities
├── python/gelexy/         # Python package
│   ├── __init__.py        # Package initialization
│   ├── model/             # Python model classes
│   ├── data/              # Data handling utilities
│   ├── predictor/         # Prediction utilities
│   └── utils/             # Python utilities
├── tests/                 # Test suite
│   ├── python/            # Python tests
│   └── cpp/               # C++ tests
├── docs/                  # Documentation
├── notebooks/             # Jupyter notebooks
├── benchmark/             # Performance benchmarks
└── ext/                   # External dependencies
    ├── eigen/             # Eigen library
    └── mio/               # Memory-mapped I/O library
```

## Build System & Dependencies

### Core Dependencies

- **Python**: >=3.11 (3.13.3+ recommended)
- **C++**: C++20 standard
- **BLAS/LAPACK**: MKL or OpenBLAS
- **Build Tools**: CMake >=3.18, Ninja

### Python Dependencies

```python
numpy>=2.2.2,<3
pandas>=2.2.3,<3
formulaic>=1.1.1,<2
scipy>=1.15.2,<2
formulae>=0.5.4,<0.6
h5py>=3.12.1,<4
```

### C++ Dependencies

- **nanobind**: Python binding generation
- **spdlog**: Logging library
- **Eigen**: Linear algebra
- **Armadillo**: Numerical library
- **mio**: Memory-mapped I/O

### Build Configuration

The project uses **scikit-build-core** with CMake. Key build options:

- `USE_MKL=ON/OFF`: Use Intel MKL or OpenBLAS
- `BUILD_PYTHON=ON/OFF`: Build Python bindings
- `BUILD_TEST=ON/OFF`: Build tests

## Core Architecture

### C++ Core Components

1. **Model Layer** (`include/gelex/model/`, `src/model/`)
   - Bayesian model implementations
   - Effect management system
   - MCMC sampling algorithms

2. **Data Layer** (`include/gelex/data/`, `src/data/`)
   - GRM (Genomic Relationship Matrix) handling
   - BED file reader for genotype data
   - Design matrix construction

3. **Optimization Layer** (`include/gelex/optim/`, `src/optim/`)
   - Optimization policies
   - Estimator implementations

### Python Bindings (`src/python/`)

- **Bayes bindings** (`bayes.cpp`): Bayesian model interfaces
- **Frequentist bindings** (`freq.cpp`): GBLUP interfaces
- **Core bindings** (`bindings.cpp`): Main binding infrastructure

### Python Package (`python/gelexy/`)

- **Model classes**: GBLUP, BayesModel, ModelMaker
- **Data utilities**: GRM loading, phenotype reading
- **Predictors**: Bayesian prediction utilities
- **Utils**: Cross-validation, alignment, timing

## Key Models & Algorithms

### Bayesian Models Supported

- **BayesA**: Single-component Bayes A
- **BayesRR**: Ridge regression Bayes
- **BayesB**: Two-component mixture
- **BayesBpi**: Bayes B with variable π
- **BayesC**: Multi-component mixture
- **BayesCpi**: Bayes C with variable π
- **BayesR**: Multi-component with fixed variances

### GBLUP Model

- Standard genomic best linear unbiased prediction
- Support for multiple random effects
- Efficient matrix operations

### MCMC Sampling

- Configurable burn-in, thinning, and chain parameters
- Posterior sampling and summary statistics
- Convergence diagnostics (ESS, R-hat)

## Data Formats & I/O

### Supported Formats

- **BED**: PLINK binary genotype format
- **TSV**: Tab-separated phenotype data
- **HDF5**: Model parameter storage
- **GRM**: Binary genomic relationship matrices

### Data Processing

- Genotype-phenotype alignment
- Missing value handling
- Categorical variable conversion
- Design matrix construction

## Testing Framework

### Python Tests (`tests/python/`)

- Model initialization and properties
- Formula parsing
- GRM functionality
- Bayesian priors
- Alignment utilities

### C++ Tests (`tests/cpp/`)

- GRM calculations
- Bed file reading
- Design matrix construction
- MCMC parameter testing
- Mathematical utilities

### Test Data

Located in `tests/data/` with sample genotype and phenotype files.

## Development Workflow

### Environment Setup

```bash
# Using pixi (recommended)
pixi install
pixi run build

# Using conda/mamba
mamba env create -f environment.yaml
mamba activate gelexy
pip install -e .
```

### Building from Source

```bash
# Configure with MKL
cmake -B build -DUSE_MKL=ON -DBUILD_PYTHON=ON -DBUILD_TEST=ON

# Build
cmake --build build -j16

# Run tests
cd build && ctest --output-on-failure
```

### Python Development

```bash
# Install in development mode
pip install -e .

# Run Python tests
pytest tests/python/

# Build wheel
uv build
```

## Performance Considerations

### BLAS Configuration

- **MKL**: Better performance on Intel systems
- **OpenBLAS**: Good cross-platform performance
- Configure with `USE_MKL=ON/OFF` in CMake

### Memory Management

- Uses memory-mapped I/O for large genotype files
- Efficient matrix operations with Eigen/Armadillo
- Sparse matrix support for design matrices

### Parallelization

- OpenMP support for multi-core operations
- Thread-safe implementations where appropriate

## Key Files & Entry Points

### Main Python Interface

- `python/gelexy/__init__.py`: Primary package exports
- `python/gelexy/model/__init__.py`: Model factory functions
- `python/gelexy/data/__init__.py`: Data loading utilities

### Core C++ Components

- `src/python/bindings.cpp`: Main binding registration
- `src/model/effects_manager.h`: Effect management system
- `src/data/bed_reader.cpp`: Genotype file reader
- `src/data/grm.cpp`: GRM calculations

### Configuration Files

- `pyproject.toml`: Python package configuration
- `CMakeLists.txt`: C++ build configuration
- `environment.yaml`: Conda environment
- `pixi.toml`: Pixi task runner configuration

## Common Tasks & Examples

### GBLUP Example

```python
import gelexy as gx

# Create model maker
mk_m = gx.make_gblup("phenotype_data.tsv")

# Create GRMs
a = gx.make_grm("genotype_data.bed", method="additive")
d = gx.make_grm("genotype_data.bed", method="dominant")

# Build and fit model
model = mk_m.make("Yield ~ 1 + g[a] + g[d]", {"a": a, "d": d})
est = gx.Estimator(max_iter=20)
est.fit(model)
```

### Bayesian Example

```python
import gelexy as gx

# Create Bayesian model
model = gx.make_bayes("phenotype_data.tsv")

# Add effects and run MCMC
model.add_random_effect("genetic", grm_matrix)
model.set_sigma_prior("genetic", nu=4, s2=1.0)

mcmc = gx.MCMC(gx.MCMCParams(n_iters=5000, n_burnin=3000))
result = mcmc.run(model, seed=42)
```

## Troubleshooting & Common Issues

### Build Issues

1. **Missing BLAS/LAPACK**: Ensure MKL or OpenBLAS is installed
2. **Python version**: Requires Python 3.11+
3. **CMake version**: Requires CMake 3.18+

### Runtime Issues

1. **Memory**: Large datasets may require significant RAM
2. **File permissions**: Ensure read access to data files
3. **Matrix dimensions**: Verify genotype-phenotype alignment

### Performance Tips

1. Use MKL for best performance on Intel systems
2. Enable OpenMP for multi-core operations
3. Use memory-mapped I/O for large files

## Future Development Directions

### Immediate Priorities

1. **Documentation**: Complete API documentation
2. **Testing**: Expand test coverage
3. **Performance**: Optimize critical paths

### Medium-term Goals

1. **Additional Models**: Support for more Bayesian methods
2. **GPU Acceleration**: CUDA/OpenCL support
3. **Web Interface**: Basic web dashboard

### Long-term Vision

1. **Cloud Integration**: AWS/Azure deployment options
2. **Pipeline Integration**: Nextflow/Snakemake support
3. **Community Features**: Plugin system for custom models

## Contact & Support

- **Maintainer**: r1cheu (chenrulei@cemps.ac.cn)
- **Repository**: https://github.com/r1cheu/gelexy
- **Issues**: GitHub issue tracker

---

_This handover document was generated on 2025-08-26. Please keep it updated as the project evolves._
