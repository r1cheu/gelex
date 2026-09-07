// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/genetic/types.h"
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
    auto proportion = gelex::default_proportion(Modes);
    proportion.additive
        = config.genetic_variance_proportion.get<gelex::GeneticMode::A>()
              .value_or(proportion.additive);
    proportion.dominance
        = config.genetic_variance_proportion.get<gelex::GeneticMode::D>()
              .value_or(proportion.dominance);
    proportion.random = config.random_pve.value_or(0.0);
    return gelex::VarianceBudget{proportion};
}

template <typename Spec, gelex::GeneticMode Mode>
auto make_mode_spec(const McmcConfig& config) -> Spec
{
    const auto& values = config.mixture_probabilities.template get<Mode>();
    const auto probability_option
        = option_name<Mode>({.additive = "--pi", .dominance = "--dpi"});
    if constexpr (requires(const Spec& spec) { spec.probability(); })
    {
        if (values.empty())
        {
            return Spec{};
        }
        if (values.size() != 1)
        {
            throw gelex::GelexException(
                fmt::format(
                    "{} requires 1 value, got {}",
                    probability_option,
                    values.size()));
        }
        return Spec{values.front()};
    }
    else
    {
        const auto& scale_values = config.mixture_scales.template get<Mode>();
        const auto scale_option = option_name<Mode>(
            {.additive = "--scale", .dominance = "--dscale"});
        const Spec defaults;
        const auto probabilities
            = values.empty()
                  ? defaults.probabilities()
                  : to_array<Spec::class_count>(values, probability_option);
        const auto scales
            = scale_values.empty()
                  ? defaults.scales()
                  : to_array<Spec::class_count>(scale_values, scale_option);
        return Spec{probabilities, scales};
    }
}

template <gelex::GeneticModeSet Modes, gelex::BayesMethod Method>
auto make_mcmc_recipe(const McmcConfig& config)
{
    using Recipe = gelex::BuiltinBayesRecipe<Modes, Method>;
    using Spec = typename Recipe::genetic_spec_type;
    auto genetic_spec = [&]() -> Spec
    {
        if constexpr (
            Method == gelex::BayesMethod::RR || Method == gelex::BayesMethod::A)
        {
            return {};
        }
        else if constexpr (Method == gelex::BayesMethod::CD)
        {
            using JointSpec = gelex::JointSpikeSlabSpec<>;
            const auto& values = config.mixture_probabilities.joint();
            auto joint_spec = values.empty()
                                  ? JointSpec{}
                                  : JointSpec{to_array<JointSpec::class_count>(
                                        values, "--jpi")};
            return Spec{{}, std::move(joint_spec)};
        }
        else
        {
            return gelex::generate_mode_values<Modes>(
                [&]<gelex::GeneticMode Mode>()
                {
                    return make_mode_spec<
                        typename Spec::template mode_value_type<Mode>,
                        Mode>(config);
                });
        }
    }();
    return Recipe{std::move(genetic_spec), make_variance_budget<Modes>(config)};
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
            if constexpr (Modes == additive_dominance_mode)
            {
                return std::invoke(
                    std::forward<Function>(function),
                    make_mcmc_recipe<Modes, gelex::BayesMethod::CD>(config));
            }
            else
            {
                throw gelex::GelexException("--method CD requires --mode AD");
            }
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
