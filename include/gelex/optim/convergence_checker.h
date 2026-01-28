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

#ifndef GELEX_OPTIM_CONVERGENCE_CHECKER_H_
#define GELEX_OPTIM_CONVERGENCE_CHECKER_H_

#include <Eigen/Dense>

namespace gelex
{
class ConvergenceChecker
{
   public:
    explicit ConvergenceChecker(double tol = 1e-8) : tol_{tol} {}

    bool is_converged(
        const Eigen::Ref<const Eigen::VectorXd>& new_sigma,
        double new_loglike);

    void clear()
    {
        converged_ = false;
        old_sigma_.resize(0);
        old_loglike_ = 0.0;
    }

   private:
    double tol_;
    bool converged_{false};
    Eigen::VectorXd old_sigma_;
    double old_loglike_{};

    void update(const Eigen::Ref<const Eigen::VectorXd>& sigma, double loglike)
    {
        old_sigma_ = sigma;
        old_loglike_ = loglike;
    }

    double compute_sigma_diff(
        const Eigen::Ref<const Eigen::VectorXd>& new_sigma) const;
    double compute_loglike_diff(double new_loglike) const;
};
}  // namespace gelex

#endif  // GELEX_OPTIM_CONVERGENCE_CHECKER_H_
