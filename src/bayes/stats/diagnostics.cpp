// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/stats/diagnostics.h"

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
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
using Eigen::VectorXd;

namespace
{

constexpr double nan = std::numeric_limits<double>::quiet_NaN();

// Smallest number >= target whose only prime factors are 2, 3 and 5, like
// scipy.fftpack.next_fast_len.
auto fft_next_fast_len(Index target) -> Index
{
    if (target <= 2)
    {
        return target;
    }
    while (true)
    {
        Index m = target;
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

// Biased autocovariance of one series via FFT: lag k is
// sum_t (x_t - mean)(x_{t+k} - mean) / n.
auto autocovariance(const ChainsView::ConstColXpr& x) -> VectorXd
{
    const Index n_draws = x.size();
    const Index M2 = 2 * fft_next_fast_len(n_draws);

    VectorXd padding_signal = VectorXd::Zero(M2);
    padding_signal.head(n_draws) = x.array() - x.mean();

    Eigen::FFT<double> fft;
    Eigen::VectorXcd freqvec;
    fft.fwd(freqvec, padding_signal);

    Eigen::VectorXcd freqvec_gram
        = freqvec.array() * freqvec.conjugate().array();

    Eigen::VectorXcd autocov_cx;
    fft.inv(autocov_cx, freqvec_gram);

    return autocov_cx.real().head(n_draws) / static_cast<double>(n_draws);
}

// Geyer's initial monotone sequence: sum adjacent lag pairs, clip at zero,
// force monotone decrease, and return 2 * sum - 1 (the integrated
// autocorrelation time).
auto geyer_sum(const VectorXd& rho) -> double
{
    const Index n_pairs = rho.size() / 2;
    double current_min = rho(0) + rho(1);
    double sum = current_min;
    for (Index j = 1; j < n_pairs; ++j)
    {
        double val = std::max(rho(2 * j) + rho((2 * j) + 1), 0.0);
        val = std::min(val, current_min);
        current_min = val;
        sum += val;
    }
    return (2.0 * sum) - 1.0;
}

// Gelman & Rubin pooled variance over the columns of `chains`: mean
// within-chain sample variance and the between-corrected total estimate.
struct PooledVariance
{
    double within;
    double total;
};

template <typename Derived>
auto pooled_variance(const Eigen::DenseBase<Derived>& chains) -> PooledVariance
{
    const auto n_draws = static_cast<double>(chains.rows());
    const Eigen::RowVectorXd vars = matvar<0>(chains, VarNormType::Sample);
    const double within = vars.mean();
    const double population = within * (n_draws - 1.0) / n_draws;
    if (chains.cols() == 1)
    {
        return {population, population};
    }
    const Eigen::RowVectorXd means = chains.colwise().mean();
    return {within, population + vecvar(means, VarNormType::Sample)};
}

auto effective_sample_size(const ChainsView& draws) -> double
{
    const Index n_draws = draws.rows();
    const Index n_chains = draws.cols();

    VectorXd gamma_mean = VectorXd::Zero(n_draws);
    for (Index chain = 0; chain < n_chains; ++chain)
    {
        gamma_mean += autocovariance(draws.col(chain));
    }
    gamma_mean /= static_cast<double>(n_chains);

    const auto [within, total] = pooled_variance(draws);
    if (!(total > 0.0))
    {
        return nan;
    }
    VectorXd rho = 1.0 - (within - gamma_mean.array()) / total;
    rho(0) = 1.0;
    return static_cast<double>(n_chains * n_draws) / geyer_sum(rho);
}

auto split_rhat(const ChainsView& draws) -> double
{
    const Index n_half = draws.rows() / 2;
    const Index n_chains = draws.cols();
    MatrixXd halves(n_half, 2 * n_chains);
    halves.leftCols(n_chains) = draws.topRows(n_half);
    halves.rightCols(n_chains) = draws.middleRows(n_half, n_half);

    const auto [within, total] = pooled_variance(halves);
    return within > 0.0 ? std::sqrt(total / within) : nan;
}

auto hpdi_sorted(const VectorXd& sorted, double prob)
    -> std::pair<double, double>
{
    if (prob == 1)
    {
        return {sorted(0), sorted(sorted.size() - 1)};
    }
    const Index mass = sorted.size();
    const auto index_length
        = static_cast<Index>(prob * static_cast<double>(mass));
    const Index tails = mass - index_length;

    const VectorXd intervals
        = sorted.tail(tails).array() - sorted.head(tails).array();

    Index index_start{};
    intervals.minCoeff(&index_start);

    return {sorted(index_start), sorted(index_start + index_length)};
}

auto median_sorted(const VectorXd& sorted) -> double
{
    const Index n = sorted.size();
    const Index mid = n / 2;
    return (n % 2 == 0) ? (sorted(mid - 1) + sorted(mid)) / 2.0 : sorted(mid);
}

}  // namespace

auto diagnose_chain(const ChainsView& draws, double prob) -> ChainDiagnostics
{
    if (draws.rows() < 4 || draws.cols() < 1)
    {
        throw GelexException(
            fmt::format(
                "diagnose_chain requires at least 4 draws per chain and one "
                "chain, got ({}, {})",
                draws.rows(),
                draws.cols()));
    }
    if (prob <= 0.0 || prob > 1.0)
    {
        throw GelexException(
            fmt::format("hpdi probability must lie in (0, 1], got {}", prob));
    }

    VectorXd sorted(draws.size());
    for (Index chain = 0; chain < draws.cols(); ++chain)
    {
        sorted.segment(chain * draws.rows(), draws.rows()) = draws.col(chain);
    }
    std::sort(sorted.begin(), sorted.end());
    const auto [hpdi_lower, hpdi_upper] = hpdi_sorted(sorted, prob);

    const double mean = sorted.mean();
    const double sd = std::sqrt(vecvar(sorted, VarNormType::Sample));
    const double ess = effective_sample_size(draws);

    return ChainDiagnostics{
        .mean = mean,
        .sd = sd,
        .median = median_sorted(sorted),
        .hpdi_lower = hpdi_lower,
        .hpdi_upper = hpdi_upper,
        .ess = ess,
        .mcse = sd / std::sqrt(ess),
        .split_rhat = split_rhat(draws),
    };
}

}  // namespace gelex
