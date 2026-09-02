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

#ifndef GELEX_FREQ_REML_ESTIMATOR_H_
#define GELEX_FREQ_REML_ESTIMATOR_H_

#include <cstddef>

#include "gelex/freq/model.h"
#include "gelex/freq/reml/convergence_checker.h"
#include "gelex/freq/reml/progress.h"
#include "gelex/freq/reml/summary.h"

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
        bool em_init = true) -> RemlFit;

   private:
    ConvergenceChecker convergence_checker_;
    size_t max_iter_{100};
    RemlObserver observer_;
};

}  // namespace gelex

#endif  // GELEX_FREQ_REML_ESTIMATOR_H_
