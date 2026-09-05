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

#ifndef GELEX_BAYES_GENETIC_RESULT_H_
#define GELEX_BAYES_GENETIC_RESULT_H_

#include <Eigen/Core>
#include <cstddef>
#include <utility>

#include "gelex/bayes/basic_result.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/exception.h"

namespace gelex
{

class MarkerPipResult
{
   public:
    explicit MarkerPipResult(Eigen::VectorXd probabilities)
        : probabilities_{std::move(probabilities)}
    {
        if (probabilities_.size() == 0)
        {
            throw GelexException("marker PIP probabilities are empty");
        }
        if (!probabilities_.allFinite()
            || !(probabilities_.array() >= 0.0).all()
            || !(probabilities_.array() <= 1.0).all())
        {
            throw GelexException(
                "marker PIP probabilities must be finite and within [0, 1]");
        }
    }

    [[nodiscard]] auto probabilities() const noexcept -> const Eigen::VectorXd&
    {
        return probabilities_;
    }

   private:
    Eigen::VectorXd probabilities_;
};

class MarkerPveResult
{
   public:
    explicit MarkerPveResult(Eigen::VectorXd values)
        : values_{std::move(values)}
    {
        if (values_.size() == 0)
        {
            throw GelexException("marker PVE values are empty");
        }
        if (!values_.allFinite() || !(values_.array() >= 0.0).all())
        {
            throw GelexException(
                "marker PVE values must be finite and non-negative");
        }
    }

    [[nodiscard]] auto values() const noexcept -> const Eigen::VectorXd&
    {
        return values_;
    }

   private:
    Eigen::VectorXd values_;
};

template <typename PipResult>
class MarkerEffectResult
{
   public:
    MarkerEffectResult(
        VectorPosteriorResult coefficients,
        MarkerPveResult pve,
        PipResult pip)
        : coefficients_{std::move(coefficients)},
          pve_{std::move(pve)},
          pip_{std::move(pip)}
    {
    }

    [[nodiscard]] auto coefficients() const noexcept
        -> const VectorPosteriorResult&
    {
        return coefficients_;
    }

    [[nodiscard]] auto pve() const noexcept -> const MarkerPveResult&
    {
        return pve_;
    }

    [[nodiscard]] auto pip() const noexcept -> const PipResult& { return pip_; }

   private:
    VectorPosteriorResult coefficients_;
    MarkerPveResult pve_;
    PipResult pip_;
};

template <typename PipResult>
class JointMarkerEffectResult
{
   public:
    JointMarkerEffectResult(MarkerPveResult pve, PipResult pip)
        : pve_{std::move(pve)}, pip_{std::move(pip)}
    {
    }

    [[nodiscard]] auto pve() const noexcept -> const MarkerPveResult&
    {
        return pve_;
    }

    [[nodiscard]] auto pip() const noexcept -> const PipResult& { return pip_; }

   private:
    MarkerPveResult pve_;
    PipResult pip_;
};

struct HalfNormalPosteriorResult
{
    ScalarPosteriorResult variance;
    VectorPosteriorResult probit_coefficients;
};

template <std::size_t ClassCount, MixtureWeightUpdate WeightUpdate>
struct JointSpikeSlabPosteriorResult
{
    detail::weight_result_t<WeightUpdate, VectorPosteriorResult> probabilities;
    VectorPosteriorResult component_explained_variance;
};

}  // namespace gelex

#endif  // GELEX_BAYES_GENETIC_RESULT_H_
