// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/stats/diagnostics.h"

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <fmt/format.h>
#include <limits>
#include <unsupported/Eigen/FFT>
#include <utility>

#include "gelex/exception.h"
#include "gelex/infra/var.h"

namespace gelex
{

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXd;

namespace
{

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

auto validate_chains(const Chains& x, Index min_draws) -> void
{
    if (x.empty())
    {
        throw GelexException("at least 1 chain is required");
    }
    const Index n_params = x[0].rows();
    const Index n_draws = x[0].cols();
    if (n_draws < min_draws)
    {
        throw GelexException(
            fmt::format(
                "at least {} draws are required, got {}", min_draws, n_draws));
    }
    for (const auto& chain : x)
    {
        if (chain.rows() != n_params || chain.cols() != n_draws)
        {
            throw GelexException("every chain must have the same shape");
        }
    }
}

}  // namespace

std::pair<VectorXd, VectorXd> compute_chain_variance_stats(const Chains& x)
{
    const auto n_chains = static_cast<Index>(x.size());
    const auto n_draws = x[0].cols();
    const auto n_params = x[0].rows();

    MatrixXd chain_vars(n_params, n_chains);
    MatrixXd chain_means(n_params, n_chains);

    for (Index c = 0; c < n_chains; ++c)
    {
        chain_means.col(c) = x[c].rowwise().mean();
        chain_vars.col(c) = matvar<1>(x[c], VarNormType::Sample);
    }

    VectorXd var_within = chain_vars.rowwise().mean();
    VectorXd var_estimator = var_within * (n_draws - 1) / n_draws;

    if (n_chains > 1)
    {
        MatrixXd var_between = matvar<1>(chain_means, VarNormType::Sample);
        var_estimator += var_between;
    }
    else
    {
        var_within = var_estimator;
    }

    return {var_within, var_estimator};
}

Index fft_next_fast_len(Index target)
{
    if (target <= 2)
    {
        return target;
    }
    while (true)
    {
        size_t m = target;
        while (m % 2 == 0)
        {
            m /= 2;
        }
        while (m % 3 == 0)
        {
            m /= 3;
        }
        while (m % 5 == 0)
        {
            m /= 5;
        }
        if (m == 1)
        {
            return target;
        }
        ++target;
    }
}

VectorXd gelman_rubin(const Chains& samples)
{
    validate_chains(samples, 2);
    auto [var_within, var_estimator] = compute_chain_variance_stats(samples);
    VectorXd rhat
        = (var_within.array() > 0.0)
              .select((var_estimator.array() / var_within.array()).sqrt(), nan);
    return rhat;
}

VectorXd split_gelman_rubin(const Chains& samples)
{
    validate_chains(samples, 4);
    Chains new_samples;
    new_samples.reserve(samples.size() * 2);

    const Index n_half = samples[0].cols() / 2;

    for (const auto& chain : samples)
    {
        new_samples.emplace_back(chain.leftCols(n_half));
        new_samples.emplace_back(chain.rightCols(n_half));
    }

    return gelman_rubin(new_samples);
}

MatrixXd autocorrelation(const Ref<const MatrixXd>& x, bool bias)
{
    const Index n_draws = x.cols();
    const Index n_params = x.rows();
    const Index M = fft_next_fast_len(n_draws);
    const Index M2 = 2 * M;

    MatrixXd autocorr(n_params, n_draws);

#pragma omp parallel for default(none) \
    shared(x, autocorr, n_params, n_draws, M2, bias)
    for (Index i = 0; i < n_params; ++i)
    {
        Eigen::FFT<double> fft;
        // Extract and centralize the parameter time series
        VectorXd signal = x.row(i);
        signal.array() -= signal.mean();

        VectorXd padding_signal = VectorXd::Zero(M2);
        padding_signal.head(n_draws) = signal;

        Eigen::VectorXcd freqvec;
        fft.fwd(freqvec, padding_signal);

        Eigen::VectorXcd freqvec_gram
            = freqvec.array() * freqvec.conjugate().array();

        Eigen::VectorXcd autocorr_cx;
        fft.inv(autocorr_cx, freqvec_gram);

        VectorXd autocorr_param = autocorr_cx.real().head(n_draws);

        if (!bias)
        {
            autocorr_param.array()
                /= VectorXd::LinSpaced(n_draws, static_cast<double>(n_draws), 1)
                       .array();
        }

        const double variance = autocorr_param(0);
        if (variance > 0.0)
        {
            autocorr_param /= variance;
        }
        else
        {
            autocorr_param.setConstant(
                std::numeric_limits<double>::quiet_NaN());
        }
        autocorr.row(i) = autocorr_param;
    }

    return autocorr;
}

Chains autocorrelation(const Chains& x, bool bias)
{
    Chains result;
    result.reserve(x.size());
    for (const auto& chain : x)
    {
        result.emplace_back(autocorrelation(chain, bias));
    }
    return result;
}

Chains autocovariance(const Chains& x, bool bias)
{
    Chains result = autocorrelation(x, bias);
    MatrixXd x_var(x[0].rows(), x.size());

    const auto n_chains = static_cast<Index>(result.size());

    for (Index i = 0; i < n_chains; i++)
    {
        x_var.col(i) = matvar<1>(x[i], VarNormType::Population);
    }

    for (Index i = 0; i < n_chains; ++i)
    {
        for (Eigen::Index j = 0; j < result[i].cols(); ++j)
        {
            result[i].col(j).array() *= x_var.col(i).array();
        }
    }
    return result;
}

Eigen::VectorXd effect_sample_size(const Chains& x, bool bias)
{
    validate_chains(x, 2);
    const auto n_chains = static_cast<Index>(x.size());
    const Index n_params = x[0].rows();
    const Index n_draws = x[0].cols();

    Chains gamma_k_c = autocovariance(x, bias);

    // Compute mean across chains for each parameter and lag
    Eigen::MatrixXd gamma_k_c_mean = Eigen::MatrixXd::Zero(n_params, n_draws);
    for (const auto& mat : gamma_k_c)
    {
        gamma_k_c_mean += mat;
    }
    gamma_k_c_mean /= static_cast<double>(n_chains);

    auto [var_within, var_estimator] = compute_chain_variance_stats(x);

    Eigen::MatrixXd var_within_broadcast
        = var_within * Eigen::RowVectorXd::Ones(n_draws);
    Eigen::MatrixXd var_estimator_broadcast
        = var_estimator * Eigen::RowVectorXd::Ones(n_draws);

    Eigen::MatrixXd rho_k = MatrixXd::Ones(n_params, n_draws);
    rho_k -= ((var_within_broadcast - gamma_k_c_mean).array()
              / var_estimator_broadcast.array())
                 .matrix();
    rho_k.col(0).setOnes();

    const Index n_pairs = n_draws / 2;
    Eigen::MatrixXd Rho_k(n_params, n_pairs);

    for (Index j = 0; j < n_pairs; ++j)
    {
        Rho_k.col(j) = rho_k.col(2 * j) + rho_k.col((2 * j) + 1);
    }

    Eigen::MatrixXd Rho_mono = Rho_k;

    for (Index i = 0; i < n_params; ++i)
    {
        double current_min = Rho_k(i, 0);

        for (Index j = 1; j < n_pairs; ++j)
        {
            double val = std::max(Rho_k(i, j), 0.0);
            val = std::min(val, current_min);
            current_min = val;
            Rho_mono(i, j) = val;
        }
    }

    Eigen::VectorXd Rho_sum = Rho_mono.rowwise().sum();
    Eigen::VectorXd s2 = (2.0 * Rho_sum).array() - 1.0;
    auto total_samples = static_cast<double>(n_chains * n_draws);
    Eigen::VectorXd n_eff
        = (var_estimator.array() > 0.0).select(total_samples / s2.array(), nan);

    return n_eff;
}

Eigen::VectorXd monte_carlo_standard_error(const Chains& x, bool bias)
{
    const VectorXd n_eff = effect_sample_size(x, bias);
    const auto n_chains = static_cast<Index>(x.size());
    const Index n_params = x[0].rows();
    const Index n_draws = x[0].cols();
    MatrixXd pooled(n_params, n_chains * n_draws);
    for (Index c = 0; c < n_chains; ++c)
    {
        pooled.middleCols(c * n_draws, n_draws) = x[c];
    }
    const VectorXd variance = matvar<1>(pooled, VarNormType::Sample);
    return (variance.array() / n_eff.array()).sqrt();
}

std::pair<double, double> hpdi(const Ref<const VectorXd>& samples, double prob)
{
    if (samples.size() == 0)
    {
        throw GelexException("hpdi requires at least one sample");
    }
    if (!(prob > 0.0 && prob <= 1.0))
    {
        throw GelexException(
            fmt::format("hpdi probability must lie in (0, 1], got {}", prob));
    }
    VectorXd sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    if (prob == 1)
    {
        return {sorted(0), sorted(sorted.size() - 1)};
    }
    Index mass = sorted.size();
    auto index_length = static_cast<Index>(prob * static_cast<double>(mass));
    Index tails = mass - index_length;

    VectorXd intervals
        = sorted.tail(tails).array() - sorted.head(tails).array();

    Index index_start{};
    intervals.minCoeff(&index_start);

    Index index_end = index_start + index_length;

    return {sorted(index_start), sorted(index_end)};
}

auto hpdi(const Chains& chains, double prob) -> std::pair<MatrixXd, VectorXd>
{
    const auto n_chains = static_cast<Index>(chains.size());
    const Index n_params = chains[0].rows();
    const Index n_draws = chains[0].cols();
    const Index total = n_chains * n_draws;

    MatrixXd intervals(n_params, 2);
    VectorXd medians(n_params);

    for (Index p = 0; p < n_params; ++p)
    {
        VectorXd all_draws(total);
        for (Index c = 0; c < n_chains; ++c)
        {
            all_draws.segment(c * n_draws, n_draws) = chains[c].row(p);
        }
        auto [lo, hi] = hpdi(all_draws, prob);
        intervals(p, 0) = lo;
        intervals(p, 1) = hi;

        std::sort(all_draws.begin(), all_draws.end());
        Index mid = total / 2;
        medians(p) = (total % 2 == 0)
                         ? (all_draws(mid - 1) + all_draws(mid)) / 2.0
                         : all_draws(mid);
    }

    return {std::move(intervals), std::move(medians)};
}

}  // namespace gelex
