# AGENTS.md - Gelex Development Guide

Guidance for coding agents operating in this repository.
Prefer these commands/conventions over ad-hoc choices.

- C++23
- Tests: Catch2 v3 (`catch_discover_tests`)
- Status: beta (API breaking changes are acceptable)

## Build Commands

```bash
# useful command
pixi r install-release # install release binaries to PATH
pixi r install-debug # install debug binaries to PATH
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
- Use include guards (avoid `#pragma once`).
- Guard names should be uppercase and path-derived.
- Keep public API in `include/gelex/`; avoid exposing internals.
- Use exceptions from `include/gelex/exception.h`.
- Add comments only for non-obvious intent/constraints.
- Preserve existing license headers in source files.
- Put reusable fixtures in `tests/*_fixture.{h,cpp}`.
