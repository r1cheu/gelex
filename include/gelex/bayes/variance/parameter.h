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

#ifndef GELEX_BAYES_VARIANCE_PARAMETER_H_
#define GELEX_BAYES_VARIANCE_PARAMETER_H_

#include <cmath>

#include "gelex/exception.h"

namespace gelex
{

class ScaledInvChiSqPrior
{
   public:
    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    ScaledInvChiSqPrior(double degrees_of_freedom, double scale)
        : degrees_of_freedom_{degrees_of_freedom}, scale_{scale}
    {
        if (!std::isfinite(degrees_of_freedom_) || degrees_of_freedom_ <= 0.0
            || !std::isfinite(scale_) || scale_ <= 0.0)
        {
            throw GelexException(
                "scaled inverse chi-squared prior requires finite positive "
                "degrees of freedom and scale");
        }
    }

    [[nodiscard]] auto degrees_of_freedom() const noexcept -> double
    {
        return degrees_of_freedom_;
    }

    [[nodiscard]] auto scale() const noexcept -> double { return scale_; }

   private:
    double degrees_of_freedom_;
    double scale_;
};

class VarianceParameter
{
   public:
    VarianceParameter(double initial_value, ScaledInvChiSqPrior prior)
        : initial_value_{initial_value}, prior_{prior}
    {
        if (!std::isfinite(initial_value_) || initial_value_ <= 0.0)
        {
            throw GelexException(
                "variance parameter requires a finite positive initial "
                "value");
        }
    }

    [[nodiscard]] auto initial_value() const noexcept -> double
    {
        return initial_value_;
    }

    [[nodiscard]] auto prior() const noexcept -> const ScaledInvChiSqPrior&
    {
        return prior_;
    }

   private:
    double initial_value_;
    ScaledInvChiSqPrior prior_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_PARAMETER_H_
