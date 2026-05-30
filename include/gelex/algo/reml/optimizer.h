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

#ifndef GELEX_ALGO_REML_OPTIMIZER_H_
#define GELEX_ALGO_REML_OPTIMIZER_H_

#include <Eigen/Core>

#include "gelex/algo/reml/constrain.h"
#include "gelex/algo/reml/convergence_checker.h"
#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/algo/reml/variance_calculator.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logger.h"

namespace gelex::reml
{

class Optimizer
{
   public:
    explicit Optimizer(double tol = 1e-8)
        : convergence_checker_(tol), logger_(logging::get())
    {
    }

    template <typename Policy>
    auto step(
        const gelex::FreqModel& model,
        gelex::FreqState& state,
        OptimizerState& opt_state) -> bool;

    auto is_converged() const -> bool { return converged_; }
    auto reset() -> void
    {
        convergence_checker_.clear();
        converged_ = false;
    }

   private:
    ConvergenceChecker convergence_checker_;
    bool converged_{false};
    std::shared_ptr<spdlog::logger> logger_;
};

// helper to collect variance components from FreqState into a vector
auto collect_variance_components(const gelex::FreqState& state)
    -> Eigen::VectorXd;

// helper to distribute variance components back to FreqState
auto distribute_variance_components(
    gelex::FreqState& state,
    const Eigen::Ref<const Eigen::VectorXd>& sigma) -> void;

// Template implementation
template <typename Policy>
auto Optimizer::step(
    const gelex::FreqModel& model,
    gelex::FreqState& state,
    OptimizerState& opt_state) -> bool
{
    // 1. compute V matrix from current variance components
    compute_v(model, state, opt_state.V);

    // 2. compute V^{-1} and log|V|
    opt_state.logdet_v = v_inv_logdet(opt_state.V);

    // 3. compute projection matrix P and Py
    compute_proj(model, opt_state);

    // 4. apply policy to get new variance components
    Eigen::VectorXd sigma = Policy::apply(model, state, opt_state);

    // 5. constrain variance components to be positive
    constrain(sigma, opt_state.phenotype_variance());

    // 6. distribute back to state
    distribute_variance_components(state, sigma);

    // 7. check convergence
    double loglike = compute_loglike(model, opt_state);
    converged_ = convergence_checker_.is_converged(sigma, loglike);

    return converged_;
}

}  // namespace gelex::reml

#endif  // GELEX_ALGO_REML_OPTIMIZER_H_
