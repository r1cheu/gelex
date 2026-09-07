// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_STATS_SCALED_INV_CHI2_LOG_KERNEL_H_
#define GELEX_BAYES_STATS_SCALED_INV_CHI2_LOG_KERNEL_H_

#include <cassert>
#include <cstddef>

#include "gelex/bayes/stats/scaled_inv_chi2_distribution.h"

namespace gelex
{

class ScaledInvChi2LogKernel
{
   public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    ScaledInvChi2LogKernel(
        double log_variance_coefficient,
        double inverse_variance_coefficient)
        : log_variance_coefficient_{log_variance_coefficient},
          inverse_variance_coefficient_{inverse_variance_coefficient}
    {
        assert(log_variance_coefficient_ >= 0.0);
        assert(inverse_variance_coefficient_ >= 0.0);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    [[nodiscard]] auto scaled_inv_chi2_parameters() const
        -> ScaledInvChi2Distribution<>::param_type
    {
        const double degrees_of_freedom = log_variance_coefficient_ - 2.0;
        assert(degrees_of_freedom > 0.0);
        assert(inverse_variance_coefficient_ > 0.0);
        return ScaledInvChi2Distribution<>::param_type{
            degrees_of_freedom,
            inverse_variance_coefficient_ / degrees_of_freedom};
    }

    [[nodiscard]] auto operator+(const ScaledInvChi2LogKernel& other) const
        -> ScaledInvChi2LogKernel
    {
        return ScaledInvChi2LogKernel{
            log_variance_coefficient_ + other.log_variance_coefficient_,
            inverse_variance_coefficient_
                + other.inverse_variance_coefficient_};
    }

   private:
    double log_variance_coefficient_;
    double inverse_variance_coefficient_;
};

[[nodiscard]] inline auto make_scaled_inv_chi2_prior(
    double degrees_of_freedom,
    double scale) -> ScaledInvChi2LogKernel
{
    assert(degrees_of_freedom > 0.0);
    assert(scale > 0.0);
    return ScaledInvChi2LogKernel{
        degrees_of_freedom + 2.0, degrees_of_freedom * scale};
}

[[nodiscard]] inline auto make_normal_variance_likelihood(
    std::size_t count,
    double sum_squares) -> ScaledInvChi2LogKernel
{
    assert(sum_squares >= 0.0);
    return ScaledInvChi2LogKernel{static_cast<double>(count), sum_squares};
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_SCALED_INV_CHI2_LOG_KERNEL_H_
