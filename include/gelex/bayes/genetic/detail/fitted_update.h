// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_
#define GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_

#include <Eigen/Core>
#include <array>
#include <variant>

#include "gelex/bayes/genotype/operations.h"

namespace gelex::detail
{

// A negative component index denotes an absent contribution.
[[nodiscard]] inline auto make_fitted_update(
    Eigen::Ref<Eigen::MatrixXd> fitted_values,
    Eigen::Index old_component,
    Eigen::Index new_component,
    double old_coefficient,
    double new_coefficient) -> std::
    variant<std::monostate, bayes::AxpyTarget, std::array<bayes::AxpyTarget, 2>>
{
    const auto make_target = [&](Eigen::Index component, double delta)
    { return bayes::AxpyTarget{delta, fitted_values.col(component)}; };
    if (old_component == new_component)
    {
        const double delta = new_coefficient - old_coefficient;
        if (old_component < 0 || delta == 0.0)
        {
            return std::monostate{};
        }
        return make_target(old_component, delta);
    }
    const bool remove_old = old_component >= 0 && old_coefficient != 0.0;
    const bool add_new = new_component >= 0 && new_coefficient != 0.0;
    if (remove_old && add_new)
    {
        return std::array{
            make_target(old_component, -old_coefficient),
            make_target(new_component, new_coefficient)};
    }
    if (remove_old)
    {
        return make_target(old_component, -old_coefficient);
    }
    if (add_new)
    {
        return make_target(new_component, new_coefficient);
    }
    return std::monostate{};
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_FITTED_UPDATE_H_
