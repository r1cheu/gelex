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

#ifndef GELEX_DATA_EFFECT_SAMPLER_H_
#define GELEX_DATA_EFFECT_SAMPLER_H_

#include <random>
#include <span>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct EffectSizeClass
{
    double proportion;
    double variance;
};

struct CausalEffect
{
    double additive;
    double dominance;
    int add_class;
    int dom_class;
};

class EffectSampler
{
   public:
    struct Config
    {
        std::vector<EffectSizeClass> add_classes;
        std::vector<EffectSizeClass> dom_classes;
        bool has_dominance;
        int seed;
    };

    explicit EffectSampler(Config config);

    auto sample(Eigen::Index n_snps)
        -> std::unordered_map<Eigen::Index, CausalEffect>;

   private:
    auto assign_effect_classes(
        std::span<const EffectSizeClass> classes,
        Eigen::Index count) -> std::vector<int>;

    static void validate_effect_classes(
        std::span<const EffectSizeClass> classes,
        std::string_view label);

    Config config_;
    std::mt19937_64 rng_;
};

}  // namespace gelex

#endif  // GELEX_DATA_EFFECT_SAMPLER_H_
