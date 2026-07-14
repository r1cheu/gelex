/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

auto bayes_recipe_scheme_validator() -> CLI::Validator;

auto genetic_mode_set_validator() -> CLI::Validator;

auto rint_type_validator() -> CLI::Validator;
}  // namespace cli

#endif  // APPS_CLI_VALIDATORS_H_
