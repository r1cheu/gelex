---
name: split-header
description: Separate C++ header declarations from implementations. Adds header guards, moves non-trivial non-template function bodies to .cpp. Template code stays in the header.
user_invocable: true
disable-model-invocation: true
argument-hint: <header-file-path>
---

# Split Header

Separate declarations and implementations for the given header file: $ARGUMENTS

If no argument provided, ask which header file to process.

## Rules

1. Read the header file fully before making changes.
2. Add path-derived `#ifndef`/`#define`/`#endif` header guard if missing (no `#pragma once`).
3. Classify every function body in the header:

**Keep in header:**
- Simple getters/setters (1-3 lines), `constexpr`/`consteval` functions
- **All template functions** (keep in header, move to `// --- Implementation ---` section at bottom)
- Trivial defaulted constructors/destructors
- Simple operator overloads, friend functions

**Move to .cpp:**
- Non-template functions >3 lines
- Complex constructors, functions with loops/branches/error handling
- Functions pulling heavy includes

4. **Never** split templates into `*_impl.h` — clangd cannot handle the circular include pattern. Keep all template implementations in the same header, separated by a `// --- Implementation ---` comment after the class definition.
5. Map source path: `include/gelex/<path>/<name>.h` → `src/<path>/<name>.cpp`.
6. **Ask user to confirm path before creating any new file.**
7. If `.cpp` already exists, read it first and merge.
8. Remove header-only includes that moved implementations needed (be conservative).
9. Do NOT refactor, rename, or add comments — only separate.
10. After edits, run `pixi r build-debug` to verify.
