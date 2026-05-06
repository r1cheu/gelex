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

#ifndef GELEX_MODEL_BAYES_METHOD_NEW_H_
#define GELEX_MODEL_BAYES_METHOD_NEW_H_

#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include "gelex/model/bayes/method.h"
#include "gelex/model/bayes/prior_.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

enum class DominancePolicy : std::uint8_t
{
    symmetric,
    asymmetric,
};

struct GeneticSpec
{
    GeneticMode mode{};
    VarianceSpec variance;
    std::optional<CategoricalSpec> sign;
};

struct JointSpec
{
    GeneticSpec additive;
    GeneticSpec dominance;
};

struct GeneticTerm
{
    std::variant<GeneticSpec, JointSpec> spec;
    std::optional<Mixture> mixture;
};

struct BayesConfig
{
    BayesBase base{};
    GeneticMode mode = GeneticMode::A;
    DominancePolicy dominance = DominancePolicy::symmetric;
    bool estimate_pi = false;
};

struct BayesMethod
{
    BayesConfig config;
    std::vector<GeneticTerm> genetic_terms;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_METHOD_NEW_H_
