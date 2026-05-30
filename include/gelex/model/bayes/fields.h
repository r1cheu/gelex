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

#ifndef GELEX_MODEL_BAYES_FIELDS_H_
#define GELEX_MODEL_BAYES_FIELDS_H_

#include <array>
#include <utility>

#include "gelex/model/bayes/prior_parameters.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

template <typename T>
class Variance
{
   public:
    explicit Variance(T variance) : variance_(std::move(variance)) {}

    auto variance() -> T& { return variance_; }
    auto variance() const -> const T& { return variance_; }

   private:
    T variance_;
};

template <typename T>
class Proportion
{
   public:
    explicit Proportion(T proportion) : proportion_(std::move(proportion)) {}

    auto proportion() -> T& { return proportion_; }
    auto proportion() const -> const T& { return proportion_; }

   private:
    T proportion_;
};

template <typename T>
class SampledProportionField : public Proportion<T>
{
   public:
    explicit SampledProportionField(T proportion)
        : Proportion<T>(std::move(proportion))
    {
    }

    auto assignment() -> decltype(auto)
    {
        return (this->proportion().assignment);
    }
    auto assignment() const -> decltype(auto)
    {
        return (this->proportion().assignment);
    }
};

template <typename T>
class Multiplier
{
   public:
    explicit Multiplier(T multiplier) : multiplier_(std::move(multiplier)) {}

    auto multiplier() -> T& { return multiplier_; }
    auto multiplier() const -> const T& { return multiplier_; }

   private:
    T multiplier_;
};

template <typename T>
class AssignmentField
{
   public:
    explicit AssignmentField(T assignment) : assignment_(std::move(assignment))
    {
    }

    auto assignment() -> T& { return assignment_; }
    auto assignment() const -> const T& { return assignment_; }

   private:
    T assignment_;
};

template <typename T>
class ComponentField
{
   public:
    explicit ComponentField(T component) : component_(std::move(component)) {}

    auto component() -> T& { return component_; }
    auto component() const -> const T& { return component_; }

   private:
    T component_;
};

template <typename T>
class Variances
{
   public:
    explicit Variances(std::array<T, 2> variances)
        : variances_(std::move(variances))
    {
    }

    auto variance(GeneticMode mode) -> T&
    {
        return variances_[std::to_underlying(mode)];
    }
    auto variance(GeneticMode mode) const -> const T&
    {
        return variances_[std::to_underlying(mode)];
    }

   private:
    std::array<T, 2> variances_;
};

template <typename T>
class JointVariancesField
{
   public:
    explicit JointVariancesField(T variances) : variances_(std::move(variances))
    {
    }

    auto variance(GeneticMode mode) -> SharedMarkerVariance&
    {
        return variances_.variance(mode);
    }
    auto variance(GeneticMode mode) const -> const SharedMarkerVariance&
    {
        return variances_.variance(mode);
    }
    auto variances() -> T& { return variances_; }
    auto variances() const -> const T& { return variances_; }

   private:
    T variances_;
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_FIELDS_H_
