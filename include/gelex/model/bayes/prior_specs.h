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

#ifndef GELEX_MODEL_BAYES_PRIOR_SPECS_H_
#define GELEX_MODEL_BAYES_PRIOR_SPECS_H_

#include <cstdint>
#include <utility>

#include <Eigen/Core>

namespace gelex::bayes
{

class ScaledInvChiSqPrior
{
   public:
    ScaledInvChiSqPrior() = default;
    ScaledInvChiSqPrior(double degrees_of_freedom, double scale);

    auto degrees_of_freedom() const -> double { return degrees_of_freedom_; }
    auto scale() const -> double { return scale_; }

   private:
    double degrees_of_freedom_{-2};
    double scale_{0};
};

class VarianceSpec
{
   public:
    VarianceSpec(double initial_value, ScaledInvChiSqPrior prior);

    auto initial_value() const -> double { return initial_value_; }
    auto prior() const -> const ScaledInvChiSqPrior& { return prior_; }

   private:
    double initial_value_;
    ScaledInvChiSqPrior prior_;
};

enum class MarkerVarianceScope : std::uint8_t
{
    per_marker,
    per_effect,
};

class MarkerVarianceSpec
{
   public:
    MarkerVarianceSpec(MarkerVarianceScope scope, VarianceSpec variance);

    auto scope() const -> MarkerVarianceScope { return scope_; }
    auto variance() const -> const VarianceSpec& { return variance_; }
    auto marker_variance_size(Eigen::Index num_markers) const -> Eigen::Index
    {
        switch (scope_)
        {
            case MarkerVarianceScope::per_marker:
                return num_markers;
            case MarkerVarianceScope::per_effect:
                return 1;
        }
        std::unreachable();
    }

   private:
    MarkerVarianceScope scope_{};
    VarianceSpec variance_;
};

enum class ProportionUpdate : std::uint8_t
{
    fixed,
    sampled,
};

class ProportionSpec
{
   public:
    ProportionSpec(
        Eigen::VectorXd initial_value,
        Eigen::VectorXd concentration,
        ProportionUpdate update);

    auto initial_value() const -> const Eigen::VectorXd&
    {
        return initial_value_;
    }
    auto concentration() const -> const Eigen::VectorXd&
    {
        return concentration_;
    }
    auto update() const -> ProportionUpdate { return update_; }

    auto size() const -> Eigen::Index { return initial_value_.size(); }
    auto sampled() const -> bool
    {
        return update_ == ProportionUpdate::sampled;
    }

   private:
    Eigen::VectorXd initial_value_;
    Eigen::VectorXd concentration_;
    ProportionUpdate update_{ProportionUpdate::fixed};
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_SPECS_H_
