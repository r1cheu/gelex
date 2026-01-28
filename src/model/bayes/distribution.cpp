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

#include "distribution.h"
#include <random>

#include <Eigen/Core>

namespace gelex
{
namespace detail
{
using Eigen::Ref;
using Eigen::VectorXd;
using Eigen::VectorXi;

ScaledInvChiSq::ScaledInvChiSq(const ScaledInvChiSqParams& prior_params)
    : prior_(prior_params)
{
}

ScaledInvChiSq::ScaledInvChiSq(double initial_nu, double initial_s2)
    : prior_{initial_nu, initial_s2}
{
}

void ScaledInvChiSq::compute(
    double sum_of_squared_errors,
    Eigen::Index num_observations)
{
    if (num_observations <= 0)
    {
        return;
    }

    const double posterior_nu
        = prior_.nu + static_cast<double>(num_observations);
    const double posterior_s2
        = ((prior_.nu * prior_.s2) + sum_of_squared_errors) / posterior_nu;
    posterior_ = {posterior_nu, posterior_s2};
}

void ScaledInvChiSq::compute(double single_observation_squared_error)
{
    compute(single_observation_squared_error, 1);
}

double ScaledInvChiSq::operator()(std::mt19937_64& rng) const
{
    std::chi_squared_distribution<double> chisq{posterior_.nu};
    return (posterior_.nu * posterior_.s2) / chisq(rng);
}

}  // namespace detail
}  // namespace gelex
