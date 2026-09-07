// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_VARIANCE_CALIBRATION_H_
#define GELEX_BAYES_VARIANCE_CALIBRATION_H_

#include <array>

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

}  // namespace detail

class MarkerVarianceCalibrator
{
   public:
    explicit MarkerVarianceCalibrator(
        std::array<double, all_genetic_modes.size()> base_variances);

    auto calibrate(GeneticMode mode, double initial_activity) const
        -> VarianceParameter;

   private:
    std::array<double, all_genetic_modes.size()> base_variances_;
};

auto make_marker_variance_calibrator(
    const BayesModel& model,
    const VarianceBudget& budget) -> MarkerVarianceCalibrator;

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_CALIBRATION_H_
