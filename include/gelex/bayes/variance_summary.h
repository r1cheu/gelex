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

#ifndef GELEX_BAYES_VARIANCE_SUMMARY_H_
#define GELEX_BAYES_VARIANCE_SUMMARY_H_

#include <Eigen/Core>
#include <cmath>
#include <cstddef>
#include <fmt/format.h>
#include <string_view>
#include <utility>

#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/genetic_mode.h"
#include "gelex/types/mode_values.h"

namespace gelex
{

template <GeneticModeSet Modes>
class VarianceSummary;

template <typename GeneticPrior>
[[nodiscard]] auto make_variance_summary(const BayesState<GeneticPrior>& state)
    -> VarianceSummary<BayesState<GeneticPrior>::genetic_state_type::modes>;

template <GeneticModeSet Modes>
class VarianceSummary
{
   public:
    template <GeneticMode Mode>
    [[nodiscard]] constexpr auto genetic() const noexcept -> double
    {
        return genetic_.template get<Mode>();
    }

    [[nodiscard]] constexpr auto genetic_total() const noexcept -> double
    {
        return genetic_total_;
    }

    [[nodiscard]] constexpr auto random_total() const noexcept -> double
    {
        return random_total_;
    }

    [[nodiscard]] constexpr auto residual() const noexcept -> double
    {
        return residual_;
    }

    [[nodiscard]] constexpr auto phenotypic() const noexcept -> double
    {
        return genetic_total_ + random_total_ + residual_;
    }

    template <GeneticMode Mode>
    [[nodiscard]] constexpr auto heritability() const noexcept -> double
    {
        return genetic<Mode>() / phenotypic();
    }

    [[nodiscard]] constexpr auto total_heritability() const noexcept -> double
    {
        return genetic_total_ / phenotypic();
    }

   private:
    VarianceSummary(
        HomogeneousModeValues<Modes, double> genetic,
        double genetic_total,
        double random_total,
        double residual)
        : genetic_{std::move(genetic)},
          genetic_total_{genetic_total},
          random_total_{random_total},
          residual_{residual}
    {
        genetic_.for_each([]<GeneticMode /*Mode*/>(double value)
                          { require_component(value, "genetic"); });
        require_component(genetic_total_, "total genetic");
        require_component(random_total_, "random");
        require_component(residual_, "residual");
        if (phenotypic() <= 0.0)
        {
            throw GelexException(
                "variance summary: phenotypic variance must be positive");
        }
    }

    static auto require_component(double value, std::string_view name) -> void
    {
        if (!std::isfinite(value) || value < 0.0)
        {
            throw GelexException(
                fmt::format(
                    "variance summary: {} variance must be finite and "
                    "non-negative, got {}",
                    name,
                    value));
        }
    }

    template <typename GeneticPrior>
    friend auto make_variance_summary(const BayesState<GeneticPrior>& state)
        -> VarianceSummary<BayesState<GeneticPrior>::genetic_state_type::modes>;

    HomogeneousModeValues<Modes, double> genetic_;
    double genetic_total_;
    double random_total_;
    double residual_;
};

// decltype(auto) is load-bearing: a single-mode fold is that mode's own vector
// and several modes fold into a lazy expression, so deducing by value would
// materialize a length-n vector on every draw.
template <GeneticModeSet Modes, typename GeneticState>
[[nodiscard]] auto total_genetic_value(const GeneticState& genetic)
    -> decltype(auto)
{
    return [&]<std::size_t... Index>(
               std::index_sequence<Index...>) -> decltype(auto)
    {
        return (genetic_value(genetic.template get<Modes.at(Index)>()) + ...);
    }(std::make_index_sequence<Modes.size()>{});
}

template <typename GeneticPrior>
[[nodiscard]] auto make_variance_summary(const BayesState<GeneticPrior>& state)
    -> VarianceSummary<BayesState<GeneticPrior>::genetic_state_type::modes>
{
    constexpr auto modes = BayesState<GeneticPrior>::genetic_state_type::modes;
    const auto& genetic = state.genetic();

    double random_total = 0.0;
    for (const auto& block : state.random())
    {
        random_total += detail::vecvar(
            block.fitted_values, detail::VarNormType::Population);
    }

    return VarianceSummary<modes>{
        generate_mode_values<modes>(
            [&]<GeneticMode Mode>()
            {
                return detail::vecvar(
                    genetic_value(genetic.template get<Mode>()),
                    detail::VarNormType::Population);
            }),
        detail::vecvar(
            total_genetic_value<modes>(genetic),
            detail::VarNormType::Population),
        random_total,
        state.residual().variance};
}

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_SUMMARY_H_
