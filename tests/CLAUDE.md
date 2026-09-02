# Test Organization

- New tests mirror the corresponding path below `include/gelex/`.
  `include/gelex/bayes/spec.h` maps to `tests/bayes/test_spec.cpp`, and
  `include/gelex/bayes/genetic/independent_topology.h` maps to
  `tests/bayes/genetic/test_independent_topology.cpp`.
- Keep existing legacy tests in place. Move them gradually with `git mv` when
  the corresponding area is substantially changed; do not perform unrelated
  bulk reorganizations.
- Place fixtures in the nearest common test directory of their consumers.
- Do not keep duplicate legacy and mirrored tests after a migration.
