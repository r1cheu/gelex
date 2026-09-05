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

#ifndef APPS_CLI_MCMC_RECIPE_H_
#define APPS_CLI_MCMC_RECIPE_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <fmt/format.h>
#include <functional>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/builtin_method.h"
#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/spec.h"
#include "gelex/bayes/variance/budget.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"

#include "config.h"

namespace cli
{

namespace detail
{

inline constexpr auto additive_mode
    = gelex::GeneticModeSet{gelex::GeneticMode::A};
inline constexpr auto dominance_mode
    = gelex::GeneticModeSet{gelex::GeneticMode::D};
inline constexpr auto additive_dominance_mode = McmcConfig::option_modes;

struct ModeOptionNames
{
    std::string_view additive;
    std::string_view dominance;
};

template <std::size_t Size>
auto to_array(const std::vector<double>& values, std::string_view option)
    -> std::array<double, Size>
{
    if (values.size() != Size)
    {
        throw gelex::GelexException(
            fmt::format(
                "{} requires {} values, got {}", option, Size, values.size()));
    }
    std::array<double, Size> result{};
    std::ranges::copy(values, result.begin());
    return result;
}

template <gelex::GeneticMode Mode>
constexpr auto option_name(ModeOptionNames names) noexcept -> std::string_view
{
    if constexpr (Mode == gelex::GeneticMode::A)
    {
        return names.additive;
    }
    else
    {
        static_assert(Mode == gelex::GeneticMode::D);
        return names.dominance;
    }
}

template <gelex::GeneticModeSet Modes>
auto make_variance_budget(const McmcConfig& config) -> gelex::VarianceBudget
{
    auto shares = gelex::default_shares(Modes);
    shares.additive
        = config.genetic_variance_shares.get<gelex::GeneticMode::A>().value_or(
            shares.additive);
    shares.dominance
        = config.genetic_variance_shares.get<gelex::GeneticMode::D>().value_or(
            shares.dominance);
    shares.random = config.random_pve.value_or(0.0);
    return gelex::VarianceBudget{shares};
}

template <gelex::GeneticModeSet Modes, gelex::VarianceLayout Kind>
auto make_genetic_spec(
    gelex::GaussianFamily<Kind> /*family*/,
    const McmcConfig& /*config*/) -> gelex::Gaussian
{
    return {};
}

template <
    gelex::GeneticModeSet Modes,
    gelex::VarianceLayout Kind,
    gelex::MixtureWeightUpdate WeightUpdate>
auto make_genetic_spec(
    gelex::SpikeSlabFamily<Kind, WeightUpdate> /*family*/,
    const McmcConfig& config)
{
    return gelex::generate_mode_values<Modes>(
        [&]<gelex::GeneticMode Mode>()
        {
            const auto& values
                = config.mixture_probabilities.template get<Mode>();
            if (values.empty())
            {
                return gelex::SpikeSlab{};
            }
            const auto option
                = option_name<Mode>({.additive = "--pi", .dominance = "--dpi"});
            if (values.size() != 1)
            {
                throw gelex::GelexException(
                    fmt::format(
                        "{} requires 1 value, got {}", option, values.size()));
            }
            return gelex::SpikeSlab{values.front()};
        });
}

template <gelex::GeneticModeSet Modes, gelex::MixtureWeightUpdate WeightUpdate>
auto make_genetic_spec(
    gelex::ScaledMixtureFamily<WeightUpdate> /*family*/,
    const McmcConfig& config)
{
    return gelex::generate_mode_values<Modes>(
        [&]<gelex::GeneticMode Mode>()
        {
            const auto& probability_values
                = config.mixture_probabilities.template get<Mode>();
            const auto& scale_values
                = config.mixture_scales.template get<Mode>();
            const auto probability_option
                = option_name<Mode>({.additive = "--pi", .dominance = "--dpi"});
            const auto scale_option = option_name<Mode>(
                {.additive = "--scale", .dominance = "--dscale"});
            const gelex::ScaledMixture defaults;
            const auto probabilities
                = probability_values.empty()
                      ? defaults.probabilities()
                      : to_array<gelex::ScaledMixture::class_count>(
                            probability_values, probability_option);
            const auto scales
                = scale_values.empty()
                      ? defaults.scales()
                      : to_array<gelex::ScaledMixture::class_count>(
                            scale_values, scale_option);
            return gelex::ScaledMixture{probabilities, scales};
        });
}

template <gelex::GeneticModeSet Modes, gelex::MixtureWeightUpdate WeightUpdate>
    requires(Modes == additive_dominance_mode)
auto make_genetic_spec(
    gelex::JointSpikeSlabFamily<WeightUpdate> /*family*/,
    const McmcConfig& config)
{
    auto mode_specs
        = gelex::ModeValues<Modes, gelex::Gaussian, gelex::HalfNormal>{
            gelex::Gaussian{}, gelex::HalfNormal{}};

    const auto& probability_values = config.mixture_probabilities.joint();
    auto joint_spec = probability_values.empty()
                          ? gelex::JointSpikeSlab{}
                          : gelex::JointSpikeSlab{
                                to_array<gelex::JointSpikeSlab::class_count>(
                                    probability_values, "--jpi")};
    return gelex::JointModeValues{std::move(mode_specs), std::move(joint_spec)};
}

template <gelex::GeneticModeSet Modes, gelex::BayesMethod Method>
auto make_mcmc_recipe(const McmcConfig& config)
{
    using GeneticFamily = gelex::genetic_family_t<Method>;
    auto genetic_spec = make_genetic_spec<Modes>(GeneticFamily{}, config);
    return gelex::BuiltinBayesRecipe<Modes, Method>{
        std::move(genetic_spec), make_variance_budget<Modes>(config)};
}

template <gelex::GeneticModeSet Modes, typename Function>
decltype(auto) dispatch_mcmc_method(
    const McmcConfig& config,
    Function&& function)
{
    switch (config.method)
    {
        case gelex::BayesMethod::RR:
            return std::invoke(
                std::forward<Function>(function),
                make_mcmc_recipe<Modes, gelex::BayesMethod::RR>(config));
        case gelex::BayesMethod::A:
            return std::invoke(
                std::forward<Function>(function),
                make_mcmc_recipe<Modes, gelex::BayesMethod::A>(config));
        case gelex::BayesMethod::B:
            return std::invoke(
                std::forward<Function>(function),
                make_mcmc_recipe<Modes, gelex::BayesMethod::B>(config));
        case gelex::BayesMethod::C:
            return std::invoke(
                std::forward<Function>(function),
                make_mcmc_recipe<Modes, gelex::BayesMethod::C>(config));
        case gelex::BayesMethod::R:
            return std::invoke(
                std::forward<Function>(function),
                make_mcmc_recipe<Modes, gelex::BayesMethod::R>(config));
        case gelex::BayesMethod::CD:
            throw gelex::GelexException(
                "--method CD is not supported currently");
    }
    throw gelex::GelexException("unsupported Bayesian method");
}

}  // namespace detail

template <typename Function>
decltype(auto) dispatch_mcmc_recipe(
    const McmcConfig& config,
    Function&& function)
{
    validate_mcmc_config(config);

    if (config.mode == detail::additive_mode)
    {
        return detail::dispatch_mcmc_method<detail::additive_mode>(
            config, std::forward<Function>(function));
    }
    if (config.mode == detail::dominance_mode)
    {
        return detail::dispatch_mcmc_method<detail::dominance_mode>(
            config, std::forward<Function>(function));
    }
    if (config.mode == detail::additive_dominance_mode)
    {
        return detail::dispatch_mcmc_method<detail::additive_dominance_mode>(
            config, std::forward<Function>(function));
    }
    throw gelex::GelexException("unsupported genetic mode set");
}

}  // namespace cli

#endif  // APPS_CLI_MCMC_RECIPE_H_
