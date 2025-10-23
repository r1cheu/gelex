#include "diagnostics.h"

#include <cstddef>

#include <Eigen/Core>
#include <unsupported/Eigen/FFT>

#include "../../data/math_utils.h"
#include "gelex/estimator/bayes/samples.h"

namespace gelex
{

using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::Ref;
using Eigen::VectorXd;

std::pair<VectorXd, VectorXd> compute_chain_variance_stats(const Samples& x)
{
    const auto n_chains = static_cast<Index>(x.size());
    const auto n_draws = x[0].cols();
    const auto n_params = x[0].rows();

    MatrixXd chain_vars(n_params, n_chains);
    MatrixXd chain_means(n_params, n_chains);

    for (Index c = 0; c < n_chains; ++c)
    {
        chain_means.col(c) = x[c].rowwise().mean();
        chain_vars.col(c) = detail::var(x[c], 1, 1);
    }

    VectorXd var_within = chain_vars.rowwise().mean();
    VectorXd var_estimator = var_within * (n_draws - 1) / n_draws;

    if (n_chains > 1)
    {
        MatrixXd var_between = detail::var(chain_means, 1, 1);
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

MatrixXd gelman_rubin(const Samples& samples)
{
    auto [var_within, var_estimator] = compute_chain_variance_stats(samples);
    MatrixXd rhat = (var_estimator.array() / var_within.array()).sqrt();
    return rhat;
}

MatrixXd split_gelman_rubin(const Samples& samples)
{
    Samples new_samples;
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
    Eigen::FFT<double> fft;

    for (Index i = 0; i < n_params; ++i)
    {
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

        double variance = autocorr_param(0);
        autocorr_param /= variance;
        autocorr.row(i) = autocorr_param;
    }

    return autocorr;
}

Samples autocorrelation(const Samples& x, bool bias)
{
    Samples result;
    result.reserve(x.size());
    for (const auto& chain : x)
    {
        result.emplace_back(autocorrelation(chain, bias));
    }
    return result;
}

Samples autocovariance(const Samples& x, bool bias)
{
    Samples result = autocorrelation(x, bias);
    MatrixXd x_var(x[0].rows(), x.size());

    const auto n_chains = static_cast<Index>(result.size());

    for (Index i = 0; i < n_chains; i++)
    {
        x_var.col(i) = detail::var(x[i], 0, 1);
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

Eigen::VectorXd effect_sample_size(const Samples& x, bool bias)
{
    const auto n_chains = static_cast<Index>(x.size());
    if (n_chains == 0)
    {
        throw std::invalid_argument("At least 1 chain is required");
    }

    const Index n_params = x[0].rows();
    const Index n_draws = x[0].cols();

    if (n_draws < 2)
    {
        throw std::invalid_argument("At least 2 draws are required");
    }

    Samples gamma_k_c = autocovariance(x, bias);

    // Compute mean across chains for each parameter and lag
    Eigen::MatrixXd gamma_k_c_mean(n_params, n_draws);
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
    Eigen::VectorXd n_eff = total_samples / s2.array();

    return n_eff;
}

std::pair<double, double> hpdi(Ref<VectorXd> samples, double prob)
{
    std::sort(samples.begin(), samples.end());
    if (prob == 1)
    {
        return {samples(0), samples(samples.size() - 1)};
    }
    Index mass = samples.size();
    auto index_length = static_cast<Index>(prob * static_cast<double>(mass));
    Index tails = mass - index_length;

    VectorXd intervals_left = samples.head(tails);
    VectorXd intervals_right = samples.tail(tails);

    VectorXd intervals = intervals_right - intervals_left;

    Index index_start{};
    intervals.minCoeff(&index_start);

    Index index_end = index_start + index_length;

    return {samples(index_start), samples(index_end)};
}

}  // namespace gelex
