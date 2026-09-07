// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_COMPLETION_CHOICE_DESCRIPTIONS_H_
#define APPS_CLI_COMPLETION_CHOICE_DESCRIPTIONS_H_

#include <span>
#include <string_view>
#include <utility>

namespace cli
{
using ChoiceEntry = std::pair<std::string_view, std::string_view>;

// Per-member descriptions for enum choices, keyed by the CLI type-name token
// (e.g. "coding", "model"). Grouping by type token is required because choice
// tokens collide across enums ("A" is both BayesA and the additive mode). A
// regression test pins each group's tokens to the core name-map constant they
// mirror, so this hand-written table cannot drift from the domain enums.
auto choice_group(std::string_view type_token) -> std::span<const ChoiceEntry>;
}  // namespace cli

#endif  // APPS_CLI_COMPLETION_CHOICE_DESCRIPTIONS_H_
