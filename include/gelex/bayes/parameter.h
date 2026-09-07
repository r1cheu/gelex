// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_PARAMETER_H_
#define GELEX_BAYES_PARAMETER_H_

#include "gelex/bayes/stats/scaled_inv_chi2_log_kernel.h"

namespace gelex
{

template <typename T>
struct FixedParameter
{
    T initial;
};

template <typename T, typename Prior>
struct Parameter
{
    T initial;
    Prior prior;
};

using VarianceParameter = Parameter<double, ScaledInvChi2LogKernel>;

}  // namespace gelex

#endif  // GELEX_BAYES_PARAMETER_H_
