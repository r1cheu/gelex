# Gelex Development Guide

C++23 · Catch2 v3 · Beta (breaking changes OK)

## Build

```bash
pixi r build-debug          # debug build
pixi r build-release         # release build
pixi r test                  # all tests (via ctest)
pixi r test-catch "[tag]"    # tests by Catch2 tag
pre-commit run clang-format --files <changed_files>  # format
```

**Never** use `ctest` directly — always `pixi r test` or `pixi r test-catch`.

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
- When designing the API, adhere to STL conventions.
- Use Eigen::isApprox for testing instead of element-wise comparisons.
- Prefer Eigen init-list construction (`Eigen::MatrixXd{{1,2},{3,4}}`, `Eigen::VectorXd{{1,2,3}}`) over `<<` streaming or per-cell assignment.

<!-- ooo:START -->
<!-- ooo:VERSION:0.36.0 -->
# Ouroboros — Specification-First AI Development

> Before telling AI what to build, define what should be built.
> As Socrates asked 2,500 years ago — "What do you truly know?"
> Ouroboros turns that question into an evolutionary AI workflow engine.

Most AI coding fails at the input, not the output. Ouroboros fixes this by
**exposing hidden assumptions before any code is written**.

1. **Socratic Clarity** — Question until ambiguity <= 0.2
2. **Ontological Precision** — Solve the root problem, not symptoms
3. **Evolutionary Loops** — Each evaluation cycle feeds back into better specs

```
Interview -> Seed -> Execute -> Evaluate
    ^                              |
    |                              v
    +---------- Evolutionary Loop -+
```

## ooo Commands

Each command loads its agent/MCP on-demand. Details in each skill file.

| Command | Loads |
|---------|-------|
| `ooo` | — |
| `ooo interview` | `ouroboros:socratic-interviewer` |
| `ooo seed` | `ouroboros:seed-architect` |
| `ooo run` | MCP required |
| `ooo evolve` | MCP: `evolve_step` |
| `ooo evaluate` | `ouroboros:evaluator` |
| `ooo unstuck` | `ouroboros:{persona}` |
| `ooo status` | MCP: `session_status` |
| `ooo setup` | — |
| `ooo help` | — |

## Agents

Loaded on-demand — not preloaded.

**Core**: socratic-interviewer, ontologist, seed-architect, evaluator,
wonder, reflect, advocate, contrarian, judge
**Support**: hacker, simplifier, researcher, architect
<!-- ooo:END -->
