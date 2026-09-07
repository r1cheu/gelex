// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef APPS_CLI_VALIDATORS_H_
#define APPS_CLI_VALIDATORS_H_

namespace CLI
{
class Validator;
}  // namespace CLI

namespace cli
{
auto open_unit_interval() -> CLI::Validator;

auto non_negative_number() -> CLI::Validator;

auto genotype_method_validator() -> CLI::Validator;

auto bayes_method_validator() -> CLI::Validator;

auto genetic_mode_set_validator() -> CLI::Validator;

auto rint_type_validator() -> CLI::Validator;
}  // namespace cli

#endif  // APPS_CLI_VALIDATORS_H_
