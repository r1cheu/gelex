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
