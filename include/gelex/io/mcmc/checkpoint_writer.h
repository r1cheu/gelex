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

#ifndef GELEX_IO_MCMC_CHECKPOINT_WRITER_H_
#define GELEX_IO_MCMC_CHECKPOINT_WRITER_H_

#include <random>
#include <string_view>

#include "gelex/algo/infer/mcmc/state.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

class BayesState;

auto write_checkpoint(
    const mcmc::State& state,
    const std::mt19937_64& rng,
    const bayes::LegacyBayesMethod& method,
    std::string_view prefix) -> void;

auto write_checkpoint(
    const BayesState& state,
    const std::mt19937_64& rng,
    std::string_view prefix) -> void;

}  // namespace gelex

#endif  // GELEX_IO_MCMC_CHECKPOINT_WRITER_H_
