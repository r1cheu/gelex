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

#ifndef GELEX_MODEL_BAYES_PRIOR_H_
#define GELEX_MODEL_BAYES_PRIOR_H_

#include <cstddef>
#include <cstdint>
#include <optional>
#include <variant>

#include <Eigen/Core>

#include "gelex/model/bayes/constrain_vector.h"

namespace gelex::bayes
{

struct ScaledInvChiSqPrior
{
    double degrees_of_freedom{};
    double scale{};
};

struct OldDirichletPrior
{
    Eigen::VectorXi concentration;
};

struct DirichletPrior
{
    PositiveVector<double> concentration;
};

enum class MarkerVarianceScope : std::uint8_t
{
    per_marker,
    per_effect,
};

struct VarianceSpec
{
    double initial_value{};
    ScaledInvChiSqPrior prior;
};

class MarkerVarianceSpec
{
   public:
    MarkerVarianceSpec(MarkerVarianceScope scope, VarianceSpec variance);

    auto scope() const -> MarkerVarianceScope { return scope_; };
    auto variance() const -> const VarianceSpec& { return variance_; };
    auto marker_variance_size(Eigen::Index num_markers) const -> Eigen::Index
    {
        switch (scope_)
        {
            case MarkerVarianceScope::per_marker:
                return num_markers;
            case MarkerVarianceScope::per_effect:
                return 1;
        }
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
        Simplex<double> initial_value,
        DirichletPrior prior,
        ProportionUpdate update);

    auto initial_value() const -> const Simplex<double>&
    {
        return initial_value_;
    };

    auto prior() const -> const DirichletPrior& { return prior_; };
    auto update() const -> ProportionUpdate { return update_; };

    auto size() const -> std::size_t { return initial_value_.size(); };
    auto sampled() const -> bool
    {
        return update_ == ProportionUpdate::sampled;
    };

   private:
    Simplex<double> initial_value_;
    DirichletPrior prior_;
    ProportionUpdate update_{ProportionUpdate::fixed};
};

struct OldVarianceSpec
{
    MarkerVarianceScope scope{};
    double init{};
    ScaledInvChiSqPrior prior;

    static auto make(double phenotype_variance) -> OldVarianceSpec;
    static auto make(double init, MarkerVarianceScope scope) -> OldVarianceSpec;
};

struct CategoricalSpec
{
    Eigen::VectorXd init;
    OldDirichletPrior prior;
    bool estimate = false;
};

struct SpikeSlab
{
};

struct ScaledMixture
{
    Eigen::VectorXd multiplier;
};

struct JointMixture
{
};

using VarianceStrategy = std::variant<SpikeSlab, ScaledMixture, JointMixture>;

struct BayesPolicy;

struct Mixture
{
    VarianceStrategy strategy;
    CategoricalSpec proportions;

    static auto make(const BayesPolicy&, bool estimate_pi)
        -> std::optional<Mixture>;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_PRIOR_H_
