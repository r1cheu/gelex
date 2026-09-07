// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_VARIANCE_DETAIL_CALIBRATION_H_
#define GELEX_BAYES_VARIANCE_DETAIL_CALIBRATION_H_

#include "gelex/bayes/parameter.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

class BayesModel;
class VarianceBudget;

namespace detail
{

auto make_mean_calibrated_variance_parameter(double target)
    -> VarianceParameter;

class MarkerVarianceCalibrator
{
   public:
    MarkerVarianceCalibrator(
        const BayesModel& model,
        const VarianceBudget& budget) noexcept
        : model_{&model}, budget_{&budget}
    {
    }

    auto calibrate(GeneticMode mode, double initial_activity) const
        -> VarianceParameter;

   private:
    const BayesModel* model_;
    const VarianceBudget* budget_;
};

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_DETAIL_CALIBRATION_H_
