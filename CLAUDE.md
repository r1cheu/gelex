# AGENTS.md - Gelex Development Guide

Guidance for coding agents operating in this repository.
Prefer these commands/conventions over ad-hoc choices.

- C++23
- Tests: Catch2 v3 (`catch_discover_tests`)
- Status: beta (API breaking changes are acceptable)

## Build Commands

```bash
# useful command
pixi r build-release # compile release binaries
pixi r build-debug # compile debug binaries
pixi r test # run all tests
pre-commit run clang-format --files <changed_files> # format files

```

## Code Style Guidelines

### Naming

- Classes/structs/enums: `PascalCase`
- Functions/variables/files: `snake_case`
- Private members: trailing underscore (`member_`)
- Internal constants: prefer `kPrefixName`

### Function Signatures and Types

- Use trailing return types in declarations/definitions (`auto f() -> int`).
- Prefer `std::span` and `std::string_view` for non-owning inputs.
- Prefer `Eigen::Ref<T>` / `const Eigen::Ref<const T>&` for Eigen views.
- Use `Eigen::Index` for Eigen row/column indexing.

### Eigen Function Parameters

- **Read-only**: templatize on `MatrixBase<Derived>` / `ArrayBase<Derived>` to accept expressions without copies.
- **Writable output**: take `const MatrixBase<Derived>&` and `const_cast` internally to prevent binding to temporaries; call `.derived() = ...` to assign.
- **Resizable output**: call `.derived().resize()`, not `.resize()` on the base.
- **Multiple args**: use a separate template parameter per argument to support mixed types.
- **Base class choice**: `MatrixBase` (dense matrices), `ArrayBase` (arrays), `DenseBase` (both), `EigenBase` (any, incl. sparse/permutation).
- **Non-template alternative**: `Ref<MatrixXd>` for a single fixed layout; `Ref<Matrix<Scalar,-1,-1,RowMajor>>` for row-major — avoids copies when layout matches, creates a temporary otherwise.
- Use include guards (avoid `#pragma once`).
- Guard names should be uppercase and path-derived.
- Keep public API in `include/gelex/`; avoid exposing internals.
- Use exceptions from `include/gelex/exception.h`.
- Add comments only for non-obvious intent/constraints.
- Preserve existing license headers in source files.
- Put reusable fixtures in `tests/*_fixture.{h,cpp}`.
