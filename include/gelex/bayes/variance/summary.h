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

#include <cmath>
#include <fmt/format.h>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/bayes/genetic/state.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/state.h"
#include "gelex/exception.h"
#include "gelex/genetic_mode.h"
#include "gelex/infra/var.h"

namespace gelex
{

template <GeneticModeSet Modes>
class VarianceSummary;

template <typename GeneticPrior>
[[nodiscard]] auto make_variance_summary(const BayesState<GeneticPrior>& state)
    -> VarianceSummary<GeneticPrior::modes>;

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

    [[nodiscard]] auto random() const noexcept -> std::span<const double>
    {
        return random_;
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
        std::vector<double> random,
        double residual)
        : genetic_{std::move(genetic)},
          random_{std::move(random)},
          residual_{residual}
    {
        genetic_.for_each(
            [&]<GeneticMode /*Mode*/>(double value)
            {
                require_component(value, "genetic");
                genetic_total_ += value;
            });
        for (const double value : random_)
        {
            require_component(value, "random");
            random_total_ += value;
        }
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
        -> VarianceSummary<GeneticPrior::modes>;

    HomogeneousModeValues<Modes, double> genetic_;
    std::vector<double> random_;
    double genetic_total_{0.0};
    double random_total_{0.0};
    double residual_;
};

template <typename GeneticPrior>
[[nodiscard]] auto make_variance_summary(const BayesState<GeneticPrior>& state)
    -> VarianceSummary<GeneticPrior::modes>
{
    constexpr auto modes = GeneticPrior::modes;
    const auto& genetic = state.genetic();
    std::vector<double> random;
    random.reserve(state.random().size());
    for (const auto& block : state.random())
    {
        random.push_back(vecvar(block.fitted_values, VarNormType::Population));
    }

    return VarianceSummary<modes>{
        generate_mode_values<modes>(
            [&]<GeneticMode Mode>()
            {
                return vecvar(
                    genetic_value(genetic.template get<Mode>()),
                    VarNormType::Population);
            }),
        std::move(random),
        state.residual().variance};
}

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_SUMMARY_H_
