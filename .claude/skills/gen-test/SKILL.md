---
name: gen-test
description: Generate Catch2 unit tests for a C++ source file following project conventions. Use when the user wants to create tests for a specific header or source file.
disable-model-invocation: true
---

Generate Catch2 unit tests for: $ARGUMENTS

## Project Test Conventions

- **Framework**: Catch2 v3 (`catch2/catch_test_macros.hpp`)
- **C++ standard**: C++23
- **Style**: Chromium base, Allman braces, 80 columns, 4-space indent
- **Test file location**: `tests/test_<module_name>.cpp`
- **License header**: Apache 2.0, Copyright 2026 RuLei Chen (see examples)

## Naming Patterns

- `TEST_CASE("ClassName - description", "[tag1][tag2]")`
- Tags match the module/class being tested, e.g. `[grm]`, `[effect_sampler]`
- Use `SECTION("description")` for sub-cases within a TEST_CASE

## Assertion Patterns

- `REQUIRE(expr)` for boolean checks
- `REQUIRE_NOTHROW(expr)` for happy paths
- `REQUIRE_THROWS_AS(expr, ExceptionType)` for error cases
- `REQUIRE_THAT(value, WithinAbs(expected, tolerance))` for floating point
- Import matchers: `using Catch::Matchers::WithinAbs;`

## Common Fixtures

- `BedFixture` (`bed_fixture.h`): Creates temporary PLINK BED/BIM/FAM files for genotype tests
- `FileFixture` (`file_fixture.h`): Creates temporary files and directories

## Include Patterns

- Public headers: `#include "gelex/<module>/<header>.h"`
- Test fixtures: `#include "<fixture>.h"`
- Always add `using namespace gelex;  // NOLINT`

## Workflow

1. Read the source file specified by the user
2. Identify the class/functions to test and their public API
3. Read existing tests in `tests/` for similar modules to match style
4. Generate a test file covering:
   - Construction / validation (valid and invalid inputs)
   - Core functionality (happy paths)
   - Edge cases and error handling
5. Place the test file in `tests/test_<name>.cpp`

## Reference Examples

- Simple unit test: [examples/unit-test.cpp](examples/unit-test.cpp)
- Fixture-based test: [examples/fixture-test.cpp](examples/fixture-test.cpp)
