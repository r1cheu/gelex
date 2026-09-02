# Gelex Development Guide

Beta (breaking changes OK)

## Build

```bash
pixi r init                  # configure the full tree -> project-wide compile_commands.json (clangd)
pixi r build [preset]        # build the selected CMake preset (default debug)
pixi r test                  # all core C++ tests
pixi r test-catch "[tag]"    # focused core tests; build debug first
pixi r test-ci               # CI configuration: core + CLI tests
pixi r test-install          # installed core CMake consumer
pixi r build-cli [mode]      # core + CLI (apps/) in build/cli-<mode> (default debug)
pixi r run [mode]            # build-cli then run the gelex binary (default debug)
pixi r build-python          # core + nanobind bindings (bindings/) in build/python
pixi r test-python           # build-python then pytest the bindings
pixi r coverage              # coverage report
pixi r benchmark             # run benchmarks
```

- **Never** use `ctest` directly — use the appropriate `pixi r test*` task.
- `pixi r test-catch` runs the existing debug test binary; run `pixi r build debug` first.
- Formatting hooks apply required fixes automatically; don't hand-format or run formatters directly.
- `-D` flags belong in `CMakePresets.json` (the source of truth), not in pixi tasks; the preset decides which components build. `apps/` (CLI) and `bindings/` (Python) are opt-in via `-DGELEX_BUILD_CLI`/`GELEX_BUILD_PYTHON` (default OFF).

## Naming

- Types (class/struct/enum): `PascalCase`
- Functions/variables/files: `snake_case`
- Private instance members: trailing underscore (`member_`)
- `const`/`constexpr` objects, including static members: `snake_case`
- Macros: `UPPER_SNAKE_CASE`

## Design

- Each abstraction layer must encapsulate exactly one independent axis of change. If every implementation of a layer looks the same, collapse that layer; if one layer has to carry two axes of change, split it.
- Use `class` when a type independently owns runtime invariants: keep its representation private and enforce validity through construction and controlled mutation. Use `struct` only for transparent aggregates, compile-time traits, or internal algorithm state whose invariants are explicitly owned elsewhere.
- Templates and concepts constrain types and compile-time values; they do not replace validation of runtime numeric invariants.

## Code Style

- Only comment when the logic is not obvious; names and signatures should carry the intent
- Use `noexcept` only when the function contract is unconditionally non-throwing
- No `using namespace` in headers. Avoid namespace-scope aliases in public headers unless they are part of the API.
- Trailing return types: `auto f() -> int`
- Non-owning inputs: `std::span`, `std::string_view`
- Eigen views: `Eigen::Ref<T>` / `const Eigen::Ref<const T>&`; index with `Eigen::Index`
- Explicit `#include` only — no transitive includes
- Include guards (not `#pragma once`), uppercase path-derived names
- Do not introduce extra namespaces without a clear architectural need
- Public API in `include/gelex/`; avoid exposing internals
- Throw `gelex::GelexException` (from `include/gelex/exception.h`); don't add exception subclasses or throw raw `<stdexcept>` types
- Simple getters/setters stay in headers
- Prefer constructor member initializer lists.
- Prefer `std::views::enumerate` over manual index loops when both index and value are needed
- Every file carries the Apache-2.0 header (`Copyright <year> RuLei Chen`); add it to new files, preserve on existing
- Test placement and fixture ownership follow `tests/CLAUDE.md`; Catch2 tags are domain-scoped (`[data][bed_source]`)
- When designing the API, adhere to STL conventions.
- Use Eigen::isApprox for testing instead of element-wise comparisons.
- Prefer Eigen init-list construction (`Eigen::MatrixXd{{1,2},{3,4}}`, `Eigen::VectorXd{{1,2,3}}`) over `<<` streaming or per-cell assignment.
