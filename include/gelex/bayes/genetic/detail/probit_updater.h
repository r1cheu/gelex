// Copyright 2026 RuLei Chen
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
        Eigen::Ref<Eigen::Vector2d> probit_coefficients,
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
        probit_coefficients
            = posterior_mean
              + precision_factor.matrixU().solve(standard_normal);
    }

   private:
    MultiQuadraticLogKernel prior_;
    std::normal_distribution<double> standard_normal_distribution_;
};

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_PROBIT_UPDATER_H_
