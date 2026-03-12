# AGENTS.md - Gelex Development Guide

Guidance for coding agents operating in this repository.
Prefer these commands/conventions over ad-hoc choices.

- C++23
- Tests: Catch2 v3 (`catch_discover_tests`)
- Status: beta (API breaking changes are acceptable)

## Design Principles

Prefer simple, direct solutions. Do not over-engineer. When refactoring, aim to REMOVE complexity, not add indirection.

## Code Changes

When fixing compilation or test errors, update the CALLERS to match the existing API — do NOT modify headers or add new overloads/shims unless explicitly asked.

## Git Workflow

Use `git mv` instead of regular file move operations when renaming or relocating files in the repository.

Keep commit messages concise — one short line summary. Do not write verbose multi-paragraph commit messages unless asked.

## Workflow

When implementing a user-provided plan, ask clarifying questions BEFORE starting if anything is ambiguous.

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
- Always `#include` the headers you directly use; never rely on transitive (indirect) includes.
- Use include guards (avoid `#pragma once`).
- Guard names should be uppercase and path-derived.
- Keep public API in `include/gelex/`; avoid exposing internals.
- Use exceptions from `include/gelex/exception.h`.
- Add comments only for non-obvious intent/constraints.
- Preserve existing license headers in source files.
- Put reusable fixtures in `tests/*_fixture.{h,cpp}`.

## Design Preferences

- **DO** update all callers when changing a function signature. Breaking changes are acceptable (beta status).
- **DON'T** add overloads, shims, or wrapper functions for backwards compatibility. One function, one signature.
- **DO** use `git mv` when renaming or moving files so git tracks history.
- **DO** use direct `switch`/`if` on enums or integer tags for type dispatch.
- **DON'T** introduce `std::visit` or visitor patterns unless the variant is truly open-ended.
- **DO** keep commit messages to one concise line (subject only, no body).
- **DO** fix test call sites to match the current API when tests break after refactoring.
- **DON'T** modify header APIs or add overloads just to make old test code compile.
