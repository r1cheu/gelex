// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_DETAIL_PROBIT_UPDATER_H_
#define GELEX_BAYES_GENETIC_DETAIL_PROBIT_UPDATER_H_

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <cassert>
#include <random>
#include <utility>

#include "gelex/bayes/stats/multi_quadratic_log_kernel.h"

namespace gelex::detail
{

class ProbitUpdater
{
   public:
    explicit ProbitUpdater(MultiQuadraticLogKernel prior)
        : prior_{std::move(prior)}
    {
    }

    auto update(
        Eigen::Ref<Eigen::Vector2d> annotation_coefficients,
        const MultiQuadraticLogKernel& likelihood,
        std::mt19937_64& rng) -> void
    {
        const auto posterior = prior_ + likelihood;
        standard_normal_distribution_.reset();
        const Eigen::LLT<Eigen::Matrix2d> precision_factor{
            posterior.quadratic()};
        assert(precision_factor.info() == Eigen::Success);
        const Eigen::Vector2d posterior_mean
            = precision_factor.solve(posterior.linear());
        const Eigen::Vector2d standard_normal{
            {standard_normal_distribution_(rng),
             standard_normal_distribution_(rng)}};
        annotation_coefficients
            = posterior_mean
              + precision_factor.matrixU().solve(standard_normal);
    }

   private:
    MultiQuadraticLogKernel prior_;
    std::normal_distribution<double> standard_normal_distribution_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_PROBIT_UPDATER_H_
