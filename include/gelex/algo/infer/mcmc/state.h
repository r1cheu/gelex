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

#ifndef GELEX_ALGO_INFER_MCMC_STATE_H_
#define GELEX_ALGO_INFER_MCMC_STATE_H_

#include <algorithm>
#include <cstdint>
#include <optional>
#include <ranges>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/legacy_method.h"
#include "gelex/model/bayes/legacy_prior.h"
#include "gelex/model/bayes/model.h"
#include "gelex/model/bayes/state.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

using TrackerVector = Eigen::VectorX<int8_t>;

struct Assignment
{
    Assignment(Eigen::Index num_markers, const Eigen::VectorXd& init_proportion)
        : tracker(TrackerVector::Zero(num_markers)),
          proportion(init_proportion),
          count(Eigen::VectorXi::Zero(init_proportion.size()))
    {
        count(0) = static_cast<int>(num_markers);
    }

    Assignment(
        TrackerVector tracker,
        Eigen::VectorXd proportion,
        Eigen::VectorXi count)
        : tracker(std::move(tracker)),
          proportion(std::move(proportion)),
          count(std::move(count))
    {
    }

    TrackerVector tracker;
    Eigen::VectorXd proportion;
    Eigen::VectorXi count;
};

struct ComponentAllocation
{
    ComponentAllocation(
        Eigen::Index num_markers,
        Eigen::Index num_samples,
        const Eigen::VectorXd& init_proportion)
        : assignment(num_markers, init_proportion),
          component_u(init_proportion.size() - 1),
          component_variance(Eigen::VectorXd::Zero(init_proportion.size() - 1))
    {
        for (auto& vec : component_u)
        {
            vec = Eigen::VectorXd::Zero(num_samples);
        }
    }

    ComponentAllocation(
        Assignment assignment,
        std::vector<Eigen::VectorXd> component_u,
        Eigen::VectorXd component_variance)
        : assignment(std::move(assignment)),
          component_u(std::move(component_u)),
          component_variance(std::move(component_variance))
    {
    }

    Assignment assignment;
    std::vector<Eigen::VectorXd> component_u;
    Eigen::VectorXd component_variance;
};

using MarkerAllocation = std::variant<Assignment, ComponentAllocation>;

struct LegacyGeneticState
{
    LegacyGeneticState(
        GeneticMode type,
        Eigen::VectorXd coeffs,
        Eigen::VectorXd u,
        double variance,
        Eigen::VectorXd marker_variance,
        std::optional<MarkerAllocation> group,
        std::optional<Assignment> sign)
        : type(type),
          coeffs(std::move(coeffs)),
          u(std::move(u)),
          variance(variance),
          marker_variance(std::move(marker_variance)),
          group(std::move(group)),
          sign(std::move(sign))
    {
    }

    LegacyGeneticState(
        const GeneticEffect& effect,
        const OldGeneticPrior& prior,
        GeneticMode mode);

    GeneticMode type;
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;

    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;

    std::optional<MarkerAllocation> group;
    std::optional<Assignment> sign;
};

}  // namespace gelex::bayes

namespace gelex
{

template <typename GeneticStateT>
class LegacyInferenceState
{
   public:
    LegacyInferenceState(
        const BayesModel& model,
        const bayes::LegacyBayesMethod& method);
    LegacyInferenceState(
        bayes::FixedState fixed,
        std::vector<bayes::RandomState> random,
        std::vector<GeneticStateT> genetics,
        bayes::ResidualState residual);

    auto fixed() -> bayes::FixedState& { return fixed_; }
    auto fixed() const -> const bayes::FixedState& { return fixed_; }
    auto random() -> std::vector<bayes::RandomState>& { return random_; }
    auto random() const -> const std::vector<bayes::RandomState>&
    {
        return random_;
    }

    auto genetics() const -> const std::vector<GeneticStateT>&
    {
        return genetics_;
    }
    auto genetics() -> std::vector<GeneticStateT>& { return genetics_; }

    auto genetic(GeneticMode type) const -> const GeneticStateT*
    {
        auto it = std::ranges::find(genetics_, type, &GeneticStateT::type);
        return it != genetics_.end() ? &*it : nullptr;
    }
    auto genetic(GeneticMode type) -> GeneticStateT*
    {
        auto it = std::ranges::find(genetics_, type, &GeneticStateT::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    auto residual() -> bayes::ResidualState& { return residual_; }
    auto residual() const -> const bayes::ResidualState& { return residual_; }

    auto compute_heritability() -> void;

   private:
    bayes::FixedState fixed_;
    std::vector<bayes::RandomState> random_;
    std::vector<GeneticStateT> genetics_;
    bayes::ResidualState residual_;
};

template <typename GeneticStateT>
LegacyInferenceState<GeneticStateT>::LegacyInferenceState(
    const BayesModel& model,
    const bayes::LegacyBayesMethod& method)
    : fixed_(model.fixed())
{
    for (const auto& prior : method.genetics)
    {
        std::visit(
            [&](const auto& spec)
            {
                using T = std::decay_t<decltype(spec)>;
                if constexpr (std::is_same_v<T, bayes::GeneticSpec>)
                {
                    const auto* eff = model.genetic(spec.mode);
                    if (eff == nullptr)
                    {
                        throw GelexException(
                            "BayesMethod ctor: missing genetic effect");
                    }
                    genetics_.emplace_back(*eff, prior, spec.mode);
                }
                else
                {
                    for (auto mode : {GeneticMode::A, GeneticMode::D})
                    {
                        const auto* eff = model.genetic(mode);
                        if (eff == nullptr)
                        {
                            throw GelexException(
                                "BayesMethod ctor: missing genetic effect for "
                                "JointSpec");
                        }
                        genetics_.emplace_back(*eff, prior, mode);
                    }
                }
            },
            prior.spec);
    }

    const auto& random_effects = model.random();
    random_.reserve(random_effects.size());
    for (std::size_t i = 0; i < random_effects.size(); ++i)
    {
        random_.emplace_back(random_effects[i], method.randoms[i].init);
    }

    residual_.y_adj = model.phenotype();
    residual_.variance = method.residual.init;
}

template <typename GeneticStateT>
LegacyInferenceState<GeneticStateT>::LegacyInferenceState(
    bayes::FixedState fixed,
    std::vector<bayes::RandomState> random,
    std::vector<GeneticStateT> genetics,
    bayes::ResidualState residual)
    : fixed_(std::move(fixed)),
      random_(std::move(random)),
      genetics_(std::move(genetics)),
      residual_(std::move(residual))
{
}

template <typename GeneticStateT>
auto LegacyInferenceState<GeneticStateT>::compute_heritability() -> void
{
    double sum_var = 0;

    for (const auto& rand : random_)
    {
        sum_var += rand.variance;
    }

    for (const auto& gen : genetics_)
    {
        sum_var += gen.variance;
    }

    sum_var += residual_.variance;

    for (auto& gen : genetics_)
    {
        gen.heritability = gen.variance / sum_var;
    }
}

namespace mcmc
{
using State = LegacyInferenceState<bayes::LegacyGeneticState>;
}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_STATE_H_
