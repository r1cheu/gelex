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

#ifndef GELEX_MODEL_BAYES_DISTRIBUTION_H_
#define GELEX_MODEL_BAYES_DISTRIBUTION_H_
#include <random>

#include <Eigen/Core>

namespace gelex
{
namespace detail
{

/**
 * @brief Sample from a Dirichlet distribution.
 *
 * @param alphas
 * @param rng
 */
inline Eigen::VectorXd dirichlet(
    const Eigen::Ref<Eigen::VectorXi>& alphas,
    std::mt19937_64& rng)
{
    Eigen::VectorXd pi = Eigen::VectorXd::Zero(alphas.size());
    for (int i = 0; i < alphas.size(); ++i)
    {
        std::gamma_distribution<double> gamma_dist(
            alphas(i) <= 0 ? 1 : alphas(i), 1.0);
        pi(i) = gamma_dist(rng);
    }
    return pi / (pi).sum();
}

struct ScaledInvChiSqParams
{
    double nu{};
    double s2{};
};

struct NormalParams
{
    double mean{};
    double var{};
};

class ScaledInvChiSq
{
   public:
    explicit ScaledInvChiSq(const ScaledInvChiSqParams& prior_params);
    ScaledInvChiSq(double initial_nu, double initial_s2);

    void compute(double sum_of_squared_errors, Eigen::Index num_observations);

    void compute(double single_observation_squared_error);

    double operator()(std::mt19937_64& rng) const;
    const ScaledInvChiSqParams& prior() { return prior_; }
    const ScaledInvChiSqParams& posterior() { return posterior_; }

   private:
    ScaledInvChiSqParams prior_;
    ScaledInvChiSqParams posterior_;
};
}  // namespace detail
}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_DISTRIBUTION_H_
