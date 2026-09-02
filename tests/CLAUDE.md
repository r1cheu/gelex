# Test Organization

- New tests mirror the corresponding path below `include/gelex/`.
- Keep existing legacy tests in place. Move them gradually with `git mv` when
  the corresponding area is substantially changed; do not perform unrelated
  bulk reorganizations.
- Place fixtures in the nearest common test directory of their consumers.
- Do not keep duplicate legacy and mirrored tests after a migration.
