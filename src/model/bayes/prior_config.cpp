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

#include "gelex/model/bayes/prior_config.h"

#include <algorithm>

#include <Eigen/Core>

namespace gelex
{

namespace
{

auto default_proportion(BayesBase base) -> Eigen::VectorXd
{
    switch (base)
    {
        case BayesBase::B:
        case BayesBase::C:
            return Eigen::VectorXd{{0.99, 0.01}};
        case BayesBase::R:
            return Eigen::VectorXd{{0.99, 0.005, 0.001, 0.001, 0.001}};
        case BayesBase::A:
        case BayesBase::RR:
            return Eigen::VectorXd{{0.0, 1.0}};
    }
    return Eigen::VectorXd{};
}

auto default_multiplier(BayesBase base) -> Eigen::VectorXd
{
    switch (base)
    {
        case BayesBase::R:
            return Eigen::VectorXd{{0.0, 0.001, 0.01, 0.1, 1.0}};
        default:
            return Eigen::VectorXd{};
    }
}

}  // namespace

PriorSetConfig::PriorSetConfig(
    BayesMethodConfig method,
    double phenotype_variance)
    : method_(method), phenotype_variance_(phenotype_variance)
{
    genetics_.push_back(
        {GeneticMode::A,
         0.5,
         default_proportion(method.base),
         default_multiplier(method.base)});

    if (method.dominance)
    {
        genetics_.push_back(
            {GeneticMode::D,
             0.2,
             default_proportion(method.base),
             default_multiplier(method.base)});
    }
}

auto PriorSetConfig::override_proportion(
    GeneticMode type,
    std::span<const double> values) -> PriorSetConfig&
{
    if (auto* g = find_genetic(type))
    {
        g->proportion = Eigen::Map<const Eigen::VectorXd>(
            values.data(), static_cast<Eigen::Index>(values.size()));
    }
    return *this;
}

auto PriorSetConfig::override_multiplier(
    GeneticMode type,
    std::span<const double> values) -> PriorSetConfig&
{
    if (auto* g = find_genetic(type))
    {
        g->multiplier = Eigen::Map<const Eigen::VectorXd>(
            values.data(), static_cast<Eigen::Index>(values.size()));
    }
    return *this;
}

auto PriorSetConfig::override_positive_prob(double value) -> PriorSetConfig&
{
    positive_prob_ = value;
    return *this;
}

auto PriorSetConfig::find_genetic(GeneticMode type) -> GeneticPriorConfig*
{
    auto it = std::ranges::find(genetics_, type, &GeneticPriorConfig::type);
    return it != genetics_.end() ? &*it : nullptr;
}

auto PriorSetConfig::find_genetic(GeneticMode type) const
    -> const GeneticPriorConfig*
{
    auto it = std::ranges::find(genetics_, type, &GeneticPriorConfig::type);
    return it != genetics_.end() ? &*it : nullptr;
}

}  // namespace gelex
