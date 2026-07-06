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

#ifndef GELEX_ALGO_REML_ESTIMATOR_H_
#define GELEX_ALGO_REML_ESTIMATOR_H_

#include <cstddef>

#include "gelex/algo/reml/optimizer.h"
#include "gelex/algo/reml/optimizer_state.h"
#include "gelex/algo/reml/result.h"
#include "gelex/freq/model.h"
#include "gelex/infra/logging/reml_event.h"

namespace gelex
{

class Estimator
{
   public:
    explicit Estimator(
        size_t max_iter = 100,
        double tol = 1e-8,
        RemlObserver observer = {});

    auto fit(
        const gelex::FreqModel& model,
        gelex::FreqState& state,
        bool em_init = true) -> RemlResult;

    auto is_converged() const -> bool { return converged_; }
    auto iter_count() const -> size_t { return iter_count_; }
    auto loglike() const -> double { return loglike_; }

   private:
    auto em_step(
        const gelex::FreqModel& model,
        gelex::FreqState& state,
        OptimizerState& opt_state) -> void;

    Optimizer optimizer_;
    size_t max_iter_{100};
    size_t iter_count_{};
    double loglike_{};
    bool converged_{false};

    RemlObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_ALGO_REML_ESTIMATOR_H_
