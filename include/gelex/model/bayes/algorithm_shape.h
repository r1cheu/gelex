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

#ifndef GELEX_MODEL_BAYES_ALGORITHM_SHAPE_H_
#define GELEX_MODEL_BAYES_ALGORITHM_SHAPE_H_

#include <cstdint>
#include <span>
#include <string_view>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

struct BayesPolicy;

enum class AlgorithmShape : std::uint8_t
{
    a_only,
    d_only,
    ad_independent,
    ad_joint,
};

[[nodiscard]] auto resolve_shape(
    const BayesPolicy& policy,
    std::span<const GeneticMode> requested) -> AlgorithmShape;

[[nodiscard]] auto to_variance_label(AlgorithmShape shape) -> std::string_view;
[[nodiscard]] auto to_heritability_label(AlgorithmShape shape)
    -> std::string_view;
[[nodiscard]] auto to_file_suffix(AlgorithmShape shape) -> std::string_view;

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_ALGORITHM_SHAPE_H_
