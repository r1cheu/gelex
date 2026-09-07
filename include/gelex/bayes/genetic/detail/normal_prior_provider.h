// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_NORMAL_PRIOR_PROVIDER_H_
#define GELEX_BAYES_GENETIC_DETAIL_NORMAL_PRIOR_PROVIDER_H_

#include <Eigen/Core>

#include "gelex/bayes/genetic/detail/marker_variance.h"
#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/stats/quadratic_log_kernel.h"

namespace gelex::detail
{

template <VarianceLayout Kind>
[[nodiscard]] auto make_normal_prior_provider(
    const marker_variance_state_t<Kind>& variance)
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return [prior = make_normal_prior(variance)](
                   Eigen::Index) -> const QuadraticLogKernel& { return prior; };
    }
    else
    {
        return [variance = &variance](Eigen::Index marker)
        { return make_normal_prior((*variance)(marker)); };
    }
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_NORMAL_PRIOR_PROVIDER_H_
