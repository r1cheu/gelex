# AGENTS.md - Gelex Development Guide

Guidance for coding agents operating in this repository.
Prefer these commands/conventions over ad-hoc choices.

## Repository Facts

- Language: C++23
- Build: CMake presets + Ninja
- Environment manager: pixi
- Tests: Catch2 v3 (`catch_discover_tests`)
- Linear algebra: Eigen (OpenBLAS default, MKL optional)
- Parallelism: OpenMP
- Status: beta (API breaking changes are acceptable)

## Build Commands

Use pixi tasks when possible (mirrors CI behavior).

```bash
# setup
pixi install

# builds
pixi run build-debug
pixi run build-release
pixi run build-native
pixi run build-reldeb
pixi run build-benchmark

# direct preset equivalents
cmake --preset debug
cmake --build --preset debug -j"$(nproc)"
cmake --preset release
cmake --build --preset release -j"$(nproc)"
```

## Test Commands

Run full suite:

```bash
pixi run test
# equivalent:
ctest --preset test-debug
```

Run a single test (preferred):

```bash
# via ctest regex
ctest --preset test-debug -R "BedPipe - load\(\) method"

# via Catch2 binary
./build/debug/tests/gelex_tests "BedPipe - load() method"
```

Discover tests/tags first:

```bash
./build/debug/tests/gelex_tests --list-tests
./build/debug/tests/gelex_tests --list-tags
./build/debug/tests/gelex_tests "[data][bed_pipe]"
ctest --preset test-debug -N
```

Coverage workflow:

```bash
pixi run config-coverage
pixi run build-coverage-bin
pixi run coverage
```

## Lint / Format / Static Analysis

Primary lint entrypoint (same as CI intent):

```bash
pre-commit run -a
```

Common targeted checks:

```bash
pre-commit run clang-format --files <changed_files>
clang-tidy -p build/debug <file.cpp>
```

If `clang-format` is not directly available in `PATH`, use the
`pre-commit` command above for formatting.

Pre-commit hooks include `clang-format`, `cmake-format`, YAML checks,
whitespace hygiene, merge-conflict checks, and debug-statement checks.

## Project Layout

- `include/gelex/`: public API headers
- `src/`: implementation + private/internal headers
- `apps/`: CLI sources (`gelex`)
- `tests/`: Catch2 suite
- `benchmark/`: benchmark targets
- `ext/`: vendored dependencies (Eigen, argparse, mio, etc.)

## Code Style Guidelines

### Formatting

- Follow `.clang-format` strictly.
- 4-space indentation, no tabs.
- Column limit: 80.
- Brace style: Allman.
- Pointer alignment: left (`Type* ptr`).

### Includes / Imports

- Include order:
  1. corresponding header (if applicable)
  2. C++ standard library
  3. External libraries (Eigen, Catch2, spdlog, OpenMP)
  4. Project headers (`gelex/...` or local internal headers)
- Keep groups separated by blank lines.

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

### Headers / API Surface

- Use include guards (avoid `#pragma once`).
- Guard names should be uppercase and path-derived.
- Keep public API in `include/gelex/`; avoid exposing internals.

### Error Handling

- Prefer exceptions from `include/gelex/exception.h`.
- Use specific exception classes when possible
  (`FileNotFoundException`, `ArgumentValidationException`, etc.).
- Include actionable context in error messages.

### Memory / Resource Management

- Use RAII.
- Prefer `std::unique_ptr` for single ownership.
- Use `std::shared_ptr` only for genuine shared lifetime.
- Avoid manual `new`/`delete`.

### Numeric / Parallel Code

- Use Eigen for linear algebra.
- Keep dimension checks explicit where safety matters.
- Use OpenMP for parallel loops; avoid data races.

### Comments and Documentation

- Add comments only for non-obvious intent/constraints.
- Avoid comments that merely restate code.
- Preserve existing license headers in source files.

### Testing Conventions

- Use Catch2 macros (`TEST_CASE`, `SECTION`, `REQUIRE`, `REQUIRE_THROWS_*`).
- Prefer descriptive test names and tags (e.g., `[data][bed_pipe]`).
- Prefer `.isApprox(...)` for Eigen numeric assertions where appropriate.
- Put reusable fixtures in `tests/*_fixture.{h,cpp}`.

## Agent Workflow Recommendations

- Prefer `pixi run ...` tasks over custom command variants.
- Run focused/single tests first, then broader suites.
- Run formatting/lint before finalizing substantial edits.
- Do not modify `CMakeLists.txt` unless explicitly requested.
