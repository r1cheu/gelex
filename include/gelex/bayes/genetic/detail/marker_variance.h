// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_MARKER_VARIANCE_H_
#define GELEX_BAYES_GENETIC_DETAIL_MARKER_VARIANCE_H_

#include <Eigen/Core>
#include <type_traits>

#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
using marker_variance_state_t = std::
    conditional_t<Kind == VarianceLayout::Pooled, double, Eigen::VectorXd>;

template <VarianceLayout Kind>
auto initial_marker_variance(
    const VarianceParameter& parameter,
    Eigen::Index marker_count) -> marker_variance_state_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return parameter.initial;
    }
    else
    {
        return Eigen::VectorXd::Constant(marker_count, parameter.initial);
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_MARKER_VARIANCE_H_
