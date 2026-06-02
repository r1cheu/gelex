# Gelex Development Guide

C++23 · Catch2 v3 · Beta (breaking changes OK)

## Build

```bash
pixi r build-debug          # debug build
pixi r build-release         # release build
pixi r test                  # all tests (via ctest)
pixi r test-catch "[tag]"    # tests by Catch2 tag
```

**Never** use `ctest` directly — always `pixi r test` or `pixi r test-catch`.

## Naming

- Types (class/struct/enum): `PascalCase`
- Functions/variables/files: `snake_case`
- Private members: trailing underscore (`member_`)
- Constants: `UPPER_SNAKE_CASE`
- Macros: `UPPER_SNAKE_CASE`

## Design

- Each abstraction layer must encapsulate exactly one independent axis of change. If every implementation of a layer looks the same, collapse that layer; if one layer has to carry two axes of change, split it.

## Code Style

- No `using namespace` in headers. Avoid namespace-scope aliases in public headers unless they are part of the API.
- Trailing return types: `auto f() -> int`
- Non-owning inputs: `std::span`, `std::string_view`
- Eigen views: `Eigen::Ref<T>` / `const Eigen::Ref<const T>&`; index with `Eigen::Index`
- Explicit `#include` only — no transitive includes
- Include order: STL, external libraries, then `gelex`; split groups with an empty line
- Prefer forward declarations when sufficient, but don't over-engineer
- Include guards (not `#pragma once`), uppercase path-derived names
- Do not introduce extra namespaces without a clear architectural need
- Public API in `include/gelex/`; avoid exposing internals
- Exceptions from `include/gelex/exception.h`
- Simple getters/setters stay in headers
- Prefer constructor member initializer lists.
- Class-private static methods over anonymous-namespace free functions
- Prefer `std::views::enumerate` over manual index loops when both index and value are needed
- Preserve existing license headers
- Reusable test fixtures in `tests/*_fixture.{h,cpp}`
- When designing the API, adhere to STL conventions.
- Use Eigen::isApprox for testing instead of element-wise comparisons.
- Prefer Eigen init-list construction (`Eigen::MatrixXd{{1,2},{3,4}}`, `Eigen::VectorXd{{1,2,3}}`) over `<<` streaming or per-cell assignment.
