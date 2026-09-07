// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

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
