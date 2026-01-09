# AGENTS.md - Gelex Development Guide

## Build & Test Commands

```bash
# Setup & build
pixi install                     # Install dependencies
pixi run build-debug            # Debug build with tests (default)
pixi run build-release          # Optimized release build
pixi run build-native           # Release with -march=native

# Testing
pixi run test                   # Run all C++ tests (ctest)
cd build/debug/tests/           # Run specific tests:
./gelex_tests "TestName*"       # Run tests matching pattern
./gelex_tests --list-tests      # List all test names
./gelex_tests [tag]             # Run tests by tag

# Installation
pixi run install-debug          # Install debug binary to ~/.local/bin
pixi run install-release        # Install release binary
```

## Code Style Guidelines

- **C++ Standard**: C++23 with CMake 3.18+
- **Naming**: PascalCase for classes, snake_case for functions/variables
- **Private members**: Use trailing underscore (variable\_)
- **Imports order**: 1) Standard library, 2) External libraries (Eigen, spdlog), 3) Internal headers
- **Error handling**: Use exceptions from `include/gelex/exception.h`
- **Memory**: RAII patterns with smart pointers (unique_ptr, shared_ptr)
- **Formatting**: Follow `.clang-format` (Chromium style, 80 columns, Allman braces)
- **Linear algebra**: Use Eigen (primary), avoid Armadillo (legacy)
- **Parallelism**: OpenMP for multi-threaded computation
- **Comments**: Only add when code is non-standard or hard to understand
- **Headers**: Public API in `include/gelex/`, implementation in `src/`
- **Data processing**: Use memory-mapped I/O (mio) with chunk-based processing

## Key Architecture Notes

- Primary linear algebra: Eigen with MKL/OpenBLAS backend
- Bayesian models: BayesAlphabet family (A, B, C, R, RR, RRD) with MCMC
- Frequentist model: GBLUP with REML
- Data handling: BED file readers, GRM computation via `DataPipe`/`BedPipe`
- Testing: Catch2 framework with sample PLINK binary files

## Coding practices

- Using span and string_view, if possible, instead of const references.
