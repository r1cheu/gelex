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

#ifndef GELEX_ALGO_SIM_EFFECT_SAMPLER_H_
#define GELEX_ALGO_SIM_EFFECT_SAMPLER_H_

#include <random>
#include <span>
#include <string_view>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct EffectSizeClass
{
    double proportion;
    double variance;
};

struct CausalEffects
{
    Eigen::VectorXd additive;
    Eigen::VectorXd dominance;
    Eigen::VectorXi add_class;
    Eigen::VectorXi dom_class;

    auto resize(Eigen::Index n_snps) -> void
    {
        additive.resize(n_snps);
        dominance.resize(n_snps);
        add_class.resize(n_snps);
        dom_class.resize(n_snps);
    }

    [[nodiscard]] auto size() const -> Eigen::Index { return additive.size(); }
};

class EffectSampler
{
   public:
    EffectSampler(
        std::vector<EffectSizeClass> add_classes,
        std::vector<EffectSizeClass> dom_classes,
        std::mt19937_64& rng);

    auto sample(Eigen::Index n_snps) -> CausalEffects;

   private:
    auto assign_effect_classes(
        std::span<const EffectSizeClass> classes,
        Eigen::Index n_snps) -> std::vector<int>;

    auto sample_effect_value(std::span<const EffectSizeClass> classes, int cls)
        -> double;

    static void validate_effect_classes(
        std::span<const EffectSizeClass> classes,
        std::string_view label);

    std::vector<EffectSizeClass> add_classes_;
    std::vector<EffectSizeClass> dom_classes_;
    std::mt19937_64& rng_;
};

}  // namespace gelex

#endif  // GELEX_ALGO_SIM_EFFECT_SAMPLER_H_
