# Gelex Development Guide

C++23 · Catch2 v3 · Beta (breaking changes OK)

## Build

```bash
pixi r build-debug          # debug build
pixi r build-release         # release build
pixi r test                  # all tests
pixi r test "[tag]"          # tests by tag
pre-commit run clang-format --files <changed_files>  # format
```

**Never** use `ctest` directly — always `pixi r test`.

## Naming

- Types (class/struct/enum): `PascalCase`
- Functions/variables/files: `snake_case`
- Private members: trailing underscore (`member_`)
- Constants: `kPrefixName`

## Code Style

- Trailing return types: `auto f() -> int`
- Non-owning inputs: `std::span`, `std::string_view`
- Eigen views: `Eigen::Ref<T>` / `const Eigen::Ref<const T>&`; index with `Eigen::Index`
- Explicit `#include` only — no transitive includes
- Include guards (not `#pragma once`), uppercase path-derived names
- Public API in `include/gelex/`; avoid exposing internals
- Exceptions from `include/gelex/exception.h`
- Simple getters/setters stay in headers
- Class-private static methods over anonymous-namespace free functions
- `switch`/`if` on enums for type dispatch
- Preserve existing license headers
- Reusable test fixtures in `tests/*_fixture.{h,cpp}`
