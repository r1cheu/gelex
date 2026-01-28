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

#include "gelex/data/grm_code_policy.h"

#include <cmath>

#include <omp.h>
#include <Eigen/Core>

namespace gelex::grm
{

namespace detail
{

struct ColStats
{
    double mean;
    Eigen::Index valid_count;
};

// Computes mean ignoring NaNs.
ColStats compute_mean_and_count(const Eigen::Ref<const Eigen::VectorXd>& col)
{
    auto is_nan = col.array().isNaN();
    Eigen::Index valid_count = col.size() - is_nan.count();

    if (valid_count == 0)
    {
        return {0.0, 0};
    }

    double sum = is_nan.select(0.0, col).sum();
    return {sum / static_cast<double>(valid_count), valid_count};
}

auto additive_mean_center(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs) -> void
{
#pragma omp parallel for default(none) shared(genotype, freqs)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        auto col = genotype.col(i);
        auto [mean, valid_count] = compute_mean_and_count(col);

        if (freqs != nullptr)
        {
            (*freqs)(i) = mean / 2.0;
        }

        if (valid_count == 0)
        {
            col.setZero();
            continue;
        }

        col = col.unaryExpr(
            [mean](double val)
            {
                if (std::isnan(val))
                {
                    return 0.0;
                }
                return val - mean;
            });
    }
}
}  // namespace detail

auto Su::operator()(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool use_additive,
    Eigen::VectorXd* freqs) const -> void
{
    if (use_additive)
    {
        detail::additive_mean_center(genotype, freqs);
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto [mean, valid_count] = detail::compute_mean_and_count(col);

            if (freqs != nullptr)
            {
                (*freqs)(i) = mean / 2.0;
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            double pA = mean / 2.0;
            double dominance_mean = 2 * pA * (1 - pA);

            col = col.unaryExpr(
                [dominance_mean](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    if (val == 2)
                    {
                        return -dominance_mean;
                    }
                    return val - dominance_mean;
                });
        }
    }
}

auto Zeng::operator()(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool use_additive,
    Eigen::VectorXd* freqs) const -> void
{
    if (use_additive)
    {
        detail::additive_mean_center(genotype, freqs);
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto [mean, valid_count] = detail::compute_mean_and_count(col);

            if (freqs != nullptr)
            {
                (*freqs)(i) = mean / 2.0;
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            double pA = mean / 2.0;
            double pa = 1 - pA;
            double v2 = -2 * pa * pa;
            double v1 = 2 * pA * pa;
            double v0 = -2 * pA * pA;

            col = col.unaryExpr(
                [v2, v1, v0](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    if (val == 2)
                    {
                        return v2;
                    }
                    if (val == 1)
                    {
                        return v1;
                    }
                    return v0;
                });
        }
    }
}

auto Yang::operator()(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool use_additive,
    Eigen::VectorXd* freqs) const -> void
{
    if (use_additive)
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto [mean, valid_count] = detail::compute_mean_and_count(col);

            if (freqs != nullptr)
            {
                (*freqs)(i) = mean / 2.0;
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            double pA = mean / 2.0;
            double pa = 1 - pA;
            double denom = std::sqrt(2 * pA * pa);

            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            double two_pa = 2 * pA;
            col = col.unaryExpr(
                [two_pa, denom](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    return (val - two_pa) / denom;
                });
        }
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto [mean, valid_count] = detail::compute_mean_and_count(col);

            if (freqs != nullptr)
            {
                (*freqs)(i) = mean / 2.0;
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            double pA = mean / 2.0;
            double pa = 1 - pA;
            double denom = 2 * pA * pa;

            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            double v2 = -2 * pa * pa / denom;
            double v1 = 2 * pA * pa / denom;
            double v0 = -2 * pA * pA / denom;

            col = col.unaryExpr(
                [v2, v1, v0](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    if (val == 2)
                    {
                        return v2;
                    }
                    if (val == 1)
                    {
                        return v1;
                    }
                    return v0;
                });
        }
    }
}

auto Vitezica::operator()(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    bool use_additive,
    Eigen::VectorXd* freqs) const -> void
{
    if (use_additive)
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto [mean, valid_count] = detail::compute_mean_and_count(col);

            if (freqs != nullptr)
            {
                (*freqs)(i) = mean / 2.0;
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            double subtract_val = mean * static_cast<double>(col.size());

            col = col.unaryExpr(
                [subtract_val](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    return val - subtract_val;
                });
        }
    }
    else
    {
#pragma omp parallel for default(none) shared(genotype, freqs)
        for (Eigen::Index i = 0; i < genotype.cols(); ++i)
        {
            auto col = genotype.col(i);
            auto is_nan = col.array().isNaN();
            Eigen::Index valid_count = col.size() - is_nan.count();

            if (freqs != nullptr)
            {
                double nAA = static_cast<double>((col.array() == 2.0).count());
                double nAa = static_cast<double>((col.array() == 1.0).count());
                (*freqs)(i) = (2.0 * nAA + nAa)
                              / (2.0 * static_cast<double>(valid_count));
            }

            if (valid_count == 0)
            {
                col.setZero();
                continue;
            }

            auto nAA = static_cast<double>((col.array() == 2.0).count());
            auto nAa = static_cast<double>((col.array() == 1.0).count());
            auto naa = static_cast<double>(valid_count) - nAA - nAa;

            double scale = static_cast<double>(col.size())
                           / static_cast<double>(valid_count);
            double f_nAA = nAA * scale;
            double f_nAa = nAa * scale;
            double f_naa = naa * scale;

            double denom = (f_nAA + f_naa - std::pow(f_nAA - f_naa, 2));
            if (denom < EPSILON)
            {
                col.setZero();
                continue;
            }

            double v2 = -2 * f_naa * f_nAa / denom;
            double v1 = 4 * f_nAA * f_naa / denom;
            double v0 = -2 * f_nAA * f_nAa / denom;

            col = col.unaryExpr(
                [v2, v1, v0](double val) -> double
                {
                    if (std::isnan(val))
                    {
                        return 0.0;
                    }
                    if (val == 2)
                    {
                        return v2;
                    }
                    if (val == 1)
                    {
                        return v1;
                    }
                    return v0;
                });
        }
    }
}

}  // namespace gelex::grm
