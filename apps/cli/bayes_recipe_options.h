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

#ifndef APPS_CLI_BAYES_RECIPE_OPTIONS_H_
#define APPS_CLI_BAYES_RECIPE_OPTIONS_H_

#include "cli/mcmc/config.h"
#include "gelex/bayes/recipe_options.h"

namespace cli
{

auto make_bayes_recipe_options(const McmcConfig& config)
    -> gelex::bayes::BayesRecipeOptions;

}  // namespace cli

#endif  // APPS_CLI_BAYES_RECIPE_OPTIONS_H_
