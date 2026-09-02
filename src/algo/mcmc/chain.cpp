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

#include "gelex/algo/mcmc/chain.h"

#include <cstddef>
#include <fmt/format.h>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/bayes/genetic/gaussian_prior.h"
#include "gelex/bayes/genetic/legacy_genetic_prior.h"
#include "gelex/bayes/legacy_prior.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{

template <typename Prior>
struct step_for_prior;

template <>
struct step_for_prior<bayes::SingleSharedGaussianPrior>
{
    using type = SingleSharedGaussianStep;
};

template <>
struct step_for_prior<bayes::SinglePerMarkerGaussianPrior>
{
    using type = SinglePerMarkerGaussianStep;
};

template <>
struct step_for_prior<bayes::SingleSharedSpikeSlabGaussianPrior>
{
    using type = SingleSharedSpikeSlabStep;
};

template <>
struct step_for_prior<bayes::SinglePerMarkerSpikeSlabGaussianPrior>
{
    using type = SinglePerMarkerSpikeSlabStep;
};

template <>
struct step_for_prior<bayes::SingleScaledMixtureGaussianPrior>
{
    using type = SingleScaledMixtureStep;
};

template <>
struct step_for_prior<bayes::JointGaussianMixturePrior>
{
    using type = JointGaussianMixtureStep;
};

template <>
struct step_for_prior<bayes::JointHalfNormalMixturePrior>
{
    using type = JointHalfNormalMixtureStep;
};

template <typename Prior>
using step_for_prior_t = typename step_for_prior<Prior>::type;

Chain::Chain(
    FixedStep fixed,
    RandomStep random,
    std::vector<SingleGeneticStep> single_genetics,
    std::vector<JointGeneticStep> joint_genetics,
    ResidualStep residual,
    BayesState& state)
    : fixed_(fixed),
      random_(random),
      single_genetics_(std::move(single_genetics)),
      joint_genetics_(std::move(joint_genetics)),
      residual_(residual),
      state_(state)
{
}

auto Chain::step() -> void
{
    fixed_.step();
    random_.step();
    for (auto& step : single_genetics_)
    {
        std::visit([](auto& value) { value.step(); }, step);
    }
    for (auto& step : joint_genetics_)
    {
        std::visit([](auto& value) { value.step(); }, step);
    }
    residual_.step();
    state_.compute_heritability();
}

auto Chain::make(
    const BayesModel& model,
    const bayes::BayesPrior& prior,
    BayesState& state,
    std::mt19937_64& rng) -> Chain
{
    auto single_genetics = std::vector<SingleGeneticStep>{};
    auto joint_genetics = std::vector<JointGeneticStep>{};
    const auto priors = prior.genetics();
    auto& blocks = state.genetics();
    if (priors.size() != blocks.size())
    {
        throw GelexException(
            fmt::format(
                "Chain::make: genetic prior/state block size mismatch: {} != "
                "{}",
                priors.size(),
                blocks.size()));
    }

    single_genetics.reserve(priors.size());
    joint_genetics.reserve(priors.size());
    for (std::size_t i = 0; i < priors.size(); ++i)
    {
        std::visit(
            [&](const auto& genetic_prior, auto& block)
            {
                using Prior = std::decay_t<decltype(genetic_prior)>;
                using Block = std::decay_t<decltype(block)>;

                if constexpr (
                    std::is_same_v<Prior, bayes::SingleGeneticPrior>
                    && std::is_same_v<Block, bayes::SingleGeneticBlockState>)
                {
                    std::visit(
                        [&]<typename Leaf>(const Leaf& leaf_prior)
                        {
                            single_genetics.emplace_back(
                                std::in_place_type<step_for_prior_t<Leaf>>,
                                model.genetic(),
                                leaf_prior,
                                block,
                                state.residual(),
                                rng);
                        },
                        genetic_prior);
                }
                else if constexpr (
                    std::is_same_v<Prior, bayes::JointGeneticPrior>
                    && std::is_same_v<Block, bayes::JointGeneticBlockState>)
                {
                    std::visit(
                        [&]<typename Leaf>(const Leaf& leaf_prior)
                        {
                            joint_genetics.emplace_back(
                                std::in_place_type<step_for_prior_t<Leaf>>,
                                model.genetic(),
                                leaf_prior,
                                block,
                                state.residual(),
                                rng);
                        },
                        genetic_prior);
                }
                else
                {
                    throw GelexException(
                        "Chain::make: genetic prior/state block type mismatch");
                }
            },
            priors[i],
            blocks[i]);
    }

    return Chain{
        FixedStep{model.fixed(), state.fixed(), state.residual(), rng},
        RandomStep{
            prior.random(),
            model.random(),
            state.random(),
            state.residual(),
            rng},
        std::move(single_genetics),
        std::move(joint_genetics),
        ResidualStep{
            model.num_individuals(), prior.residual(), state.residual(), rng},
        state};
}

}  // namespace gelex
