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

#ifndef GELEX_INFRA_STATS_DETAIL_NORMAL_H_
#define GELEX_INFRA_STATS_DETAIL_NORMAL_H_

#include <array>
#include <cmath>
#include <concepts>
#include <numbers>

#include "gelex/exception.h"

namespace gelex::detail
{
template <std::floating_point T>
auto log_norm_cdf_asymptotic(T z) -> T
{
    constexpr T LOG_2_PI = std::log(T{2} * std::numbers::pi_v<T>);
    const T z2 = z * z;
    const T z4 = z2 * z2;
    const T z6 = z4 * z2;
    const T correction = std::log1p((-T{1} / z2) + (T{3} / z4) - (T{15} / z6));
    return (-T{0.5} * z2) - (T{0.5} * LOG_2_PI) - std::log(-z) + correction;
}

template <std::floating_point T>
auto log_norm_cdf(T z) -> T
{
    constexpr T SQRT_2 = std::numbers::sqrt2_v<T>;
    const T erfc_val = std::erfc(-z / SQRT_2);
    if (erfc_val > T{0})
    {
        return std::log(T{0.5}) + std::log(erfc_val);
    }
    return log_norm_cdf_asymptotic(z);
}

inline auto norm_cdf(double x) -> double
{
    return 0.5 * (1.0 + std::erf(x / std::numbers::sqrt2));
}

inline auto norm_ppf(double p) -> double
{
    static constexpr std::array<double, 6> a
        = {-3.969683028665376e+01,
           2.209460984245205e+02,
           -2.759285104469687e+02,
           1.383577518672690e+02,
           -3.066479806614716e+01,
           2.506628277459239e+00};
    static constexpr std::array<double, 5> b
        = {-5.447609879822406e+01,
           1.615858368580409e+02,
           -1.556989798598866e+02,
           6.680131188771972e+01,
           -1.328068155288572e+01};
    static constexpr std::array<double, 6> c
        = {-7.784894002430293e-03,
           -3.223964580411365e-01,
           -2.400758277161838e+00,
           -2.549732539343734e+00,
           4.374664141464968e+00,
           2.938163982698783e+00};
    static constexpr std::array<double, 4> d
        = {7.784695709041462e-03,
           3.224671290700398e-01,
           2.445134137142996e+00,
           3.754408661907416e+00};

    static constexpr double P_LOW = 0.02425;
    static constexpr double P_HIGH = 1.0 - P_LOW;

    if (p <= 0.0 || p >= 1.0)
    {
        throw GelexException("norm_ppf: p must be in (0, 1)");
    }

    double r = 0;
    double x = 0;

    // NOLINTBEGIN(readability-math-missing-parentheses)
    if (p < P_LOW)
    {
        r = std::sqrt(-2.0 * std::log(p));
        x = (((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r
             + c[5])
            / ((((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0);
    }
    else if (p <= P_HIGH)
    {
        r = p - 0.5;
        double r2 = r * r;
        x = (((((a[0] * r2 + a[1]) * r2 + a[2]) * r2 + a[3]) * r2 + a[4]) * r2
             + a[5])
            * r
            / (((((b[0] * r2 + b[1]) * r2 + b[2]) * r2 + b[3]) * r2 + b[4]) * r2
               + 1.0);
    }
    else
    {
        r = std::sqrt(-2.0 * std::log(1.0 - p));
        x = -(((((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]) * r
              + c[5])
            / ((((d[0] * r + d[1]) * r + d[2]) * r + d[3]) * r + 1.0);
    }

    // Halley's refinement
    double e = 0.5 * std::erfc(-x * std::numbers::sqrt2_v<double> * 0.5) - p;
    double u = e * std::sqrt(2.0 * std::numbers::pi) * std::exp(0.5 * x * x);
    x -= u / (1.0 + 0.5 * x * u);
    // NOLINTEND(readability-math-missing-parentheses)

    return x;
}

}  // namespace gelex::detail

#endif  // GELEX_INFRA_STATS_DETAIL_NORMAL_H_
