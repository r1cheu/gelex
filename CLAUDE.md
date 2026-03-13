# AGENTS.md - Gelex Development Guide

Guidance for coding agents operating in this repository.
Prefer these commands/conventions over ad-hoc choices.

- C++23
- Tests: Catch2 v3 (`catch_discover_tests`)
- Status: beta (API breaking changes are acceptable)

## Build Commands

```bash
# useful command
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
- Always `#include` the headers you directly use; never rely on transitive (indirect) includes.
- Use include guards (avoid `#pragma once`).
- Guard names should be uppercase and path-derived.
- Keep public API in `include/gelex/`; avoid exposing internals.
- Use exceptions from `include/gelex/exception.h`.
- Add comments only for non-obvious intent/constraints.
- Preserve existing license headers in source files.
- Put reusable fixtures in `tests/*_fixture.{h,cpp}`.
- **DO** use `git mv` when renaming or moving files so git tracks history.
- **DO** use direct `switch`/`if` on enums or integer tags for type dispatch.
- **DO** keep commit messages to one concise line (subject only, no body).
- **DO** fix test call sites to match the current API when tests break after refactoring.
- **DON'T** modify header APIs or add overloads just to make old test code compile.
- **DO** leave simple getter setter or functions in header file.
