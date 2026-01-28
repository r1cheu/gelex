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

#ifndef GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_
#define GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_

#include <cstddef>
#include <memory>

#include "gelex/optim/optimizer.h"
#include "gelex/optim/optimizer_state.h"

namespace gelex
{

class FreqModel;
class FreqState;

namespace detail
{
class RemlLoggerBase;
}

class Estimator
{
   public:
    explicit Estimator(size_t max_iter = 100, double tol = 1e-8);
    Estimator(const Estimator&) = delete;
    Estimator(Estimator&&) = default;
    Estimator& operator=(const Estimator&) = delete;
    Estimator& operator=(Estimator&&) = default;
    Estimator(
        size_t max_iter,
        double tol,
        std::unique_ptr<detail::RemlLoggerBase> logger);
    ~Estimator();

    auto fit(
        const FreqModel& model,
        FreqState& state,
        bool em_init = true,
        bool verbose = true) -> Eigen::MatrixXd;

    auto is_converged() const -> bool { return converged_; }
    auto iter_count() const -> size_t { return iter_count_; }
    auto loglike() const -> double { return loglike_; }

   private:
    auto em_step(
        const FreqModel& model,
        FreqState& state,
        OptimizerState& opt_state) -> void;

    Optimizer optimizer_;
    size_t max_iter_{100};
    size_t iter_count_{};
    double loglike_{};
    bool converged_{false};

    std::unique_ptr<detail::RemlLoggerBase> logger_;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_ESTIMATOR_H_
