// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENETIC_PARAMETER_H_
#define GELEX_BAYES_GENETIC_PARAMETER_H_

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "gelex/bayes/genetic/types.h"
#include "gelex/bayes/parameter.h"
#include "gelex/bayes/stats/dirichlet_log_kernel.h"
#include "gelex/namespace.h"

GELEX_NAMESPACE_BEGIN(gelex)

template <MixtureWeightUpdate Update>
using ProbabilityParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<double>,
    Parameter<double, DirichletLogKernel<2>>>;

template <std::size_t ClassCount, MixtureWeightUpdate Update>
using SimplexParameter = std::conditional_t<
    Update == MixtureWeightUpdate::Disabled,
    FixedParameter<std::array<double, ClassCount>>,
    Parameter<std::array<double, ClassCount>, DirichletLogKernel<ClassCount>>>;

GELEX_NAMESPACE_BEGIN(detail)
template <MixtureWeightUpdate Update, typename T, typename Prior>
auto make_parameter(T initial, Prior prior)
{
    if constexpr (Update == MixtureWeightUpdate::Disabled)
    {
        return FixedParameter<T>{.initial = std::move(initial)};
    }
    else
    {
        return Parameter<T, Prior>{
            .initial = std::move(initial), .prior = std::move(prior)};
    }
}
GELEX_NAMESPACE_END(detail)

GELEX_NAMESPACE_END(gelex)

#endif  // GELEX_BAYES_GENETIC_PARAMETER_H_
