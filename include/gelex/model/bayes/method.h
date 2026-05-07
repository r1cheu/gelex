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

#ifndef GELEX_MODEL_BAYES_METHOD_H_
#define GELEX_MODEL_BAYES_METHOD_H_

#include <cstdint>
#include <optional>
#include <variant>
#include <vector>

#include <fmt/base.h>
#include <fmt/format.h>

#include "gelex/model/bayes/bayes_base.h"
#include "gelex/model/bayes/bayes_policy.h"
#include "gelex/model/bayes/prior.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

enum class DominancePolicy : std::uint8_t
{
    symmetric,
    asymmetric,
};

struct GeneticStats;

struct GeneticSpec
{
    GeneticMode mode{};
    VarianceSpec variance;
    std::optional<CategoricalSpec> sign;

    static auto make(const GeneticStats&, const BayesPolicy&) -> GeneticSpec;
    static auto make(const GeneticStats&, const BayesPolicy&, DominancePolicy)
        -> GeneticSpec;
};

struct JointSpec
{
    GeneticSpec additive;
    GeneticSpec dominance;
};

struct GeneticPrior
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

    constexpr auto operator==(const BayesConfig&) const -> bool = default;
};

struct BayesMethod
{
    BayesConfig config;
    std::vector<GeneticPrior> genetics;
    std::vector<VarianceSpec> randoms;
    VarianceSpec residual;
};

inline auto is_valid_method(const BayesConfig& m) -> bool
{
    if (m.base == BayesBase::kCount)
    {
        return false;
    }
    const auto& policy = policy_for(m.base);
    if (m.estimate_pi && !policy.supports_estimate_pi)
    {
        return false;
    }
    if (m.dominance == DominancePolicy::asymmetric
        && (!policy.supports_asymmetric_dominance || m.mode != GeneticMode::AD))
    {
        return false;
    }
    return true;
}

}  // namespace gelex::bayes

namespace fmt
{
template <>
struct formatter<gelex::bayes::BayesConfig> : formatter<string_view>
{
    static auto format(const gelex::bayes::BayesConfig& c, format_context& ctx)
        -> format_context::iterator
    {
        auto name = fmt::format("Bayes{}", c.base);
        if (c.estimate_pi)
        {
            name += "pi";
        }
        if (c.dominance == gelex::bayes::DominancePolicy::asymmetric)
        {
            name += " + asymmetric dominance";
        }
        else if (c.mode == gelex::GeneticMode::D)
        {
            name += " (dominance only)";
        }
        else if (c.mode == gelex::GeneticMode::AD)
        {
            name += " + dominance";
        }
        return fmt::format_to(ctx.out(), "{}", name);
    }
};

}  // namespace fmt

#endif  // GELEX_MODEL_BAYES_METHOD_H_
