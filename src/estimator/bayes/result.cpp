#include "gelex/estimator/bayes/result.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <optional>
#include <string>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "estimator/bayes/diagnostics.h"

namespace gelex
{

using Eigen::Index;
using Eigen::VectorXd;

MCMCResult::MCMCResult(MCMCSamples&& samples, double prob)
    : samples_(std::move(samples)), prob_(prob), residual_(1)
{
    fixed_
        = std::make_unique<PosteriorSummary>(samples_.fixed().coeff[0].rows());
    additive_ = samples_.additive()
                    ? std::make_unique<AdditiveSummary>(samples_.additive())
                    : nullptr;
    dominant_ = samples_.dominant()
                    ? std::make_unique<DominantSummary>(samples_.dominant())
                    : nullptr;

    for (const auto& sample : samples_.random().chain_samples)
    {
        random_.emplace_back(sample);
    }
}

void MCMCResult::compute(std::optional<double> prob)
{
    if (prob)
    {
        prob_ = prob.value();
    }

    compute_summary_statistics(*fixed_, samples_.fixed().coeff, prob_);
    for (size_t i = 0; i < samples_.random().size(); ++i)
    {
        compute_summary_statistics(
            random_[i].coeff, samples_.random().chain_samples[i].coeffs, prob_);
        compute_summary_statistics(
            random_[i].sigma, samples_.random().chain_samples[i].sigmas, prob_);
    }
    if (samples_.additive())
    {
        compute_summary_statistics(
            additive_->coeff, samples_.additive().coeffs);
        compute_summary_statistics(
            additive_->sigma, samples_.additive().sigmas);
        compute_summary_statistics(
            additive_->variance, samples_.additive().variance, prob_);
    }
    if (samples_.dominant())
    {
        compute_summary_statistics(
            dominant_->variance, samples_.dominant().variance, prob_);
    }
    compute_summary_statistics(residual_, samples_.residual(), prob_);
}

void MCMCResult::save(const std::string& prefix) const
{
    // file setting
    auto param_stream = *detail::open_file<std::ofstream>(
        prefix + ".params", std::ios_base::out);
    auto additive_stream = *detail::open_file<std::ofstream>(
        prefix + ".add", std::ios_base::out);
    std::optional<std::ofstream> dom_stream = std::nullopt;
    if (dominant_)
    {
        dom_stream = std::make_optional(*detail::open_file<std::ofstream>(
            prefix + ".dom", std::ios_base::out));
    }
    std::vector<std::string> snp_names;
    if (!samples_.bim_file_path().empty())
    {
        auto bim = *detail::BimLoader::create(samples_.bim_file_path());
        snp_names = bim.ids();
    }

    // header setting
    std::string hpdi_low
        = std::to_string(static_cast<int>(std::round(100 * (1 - prob_) / 2)));
    std::string hpdi_high
        = std::to_string(static_cast<int>(std::round(100 * (1 + prob_) / 2)));
    param_stream << "term\tmean\tstddev\t" + hpdi_low + "%" + "\t" + hpdi_high
                        + "%"
                 << "\tess\trhat\n";
    additive_stream << "snp\tmean\n";

    // function to write summary statistics
    auto write_stats = [&](const auto& stats, Index n)
    {
        for (Index i = 0; i < n; ++i)
        {
            param_stream << "\t" << stats.mean(i) << "\t" << stats.stddev(i)
                         << "\t" << stats.hpdi_low(i) << "\t"
                         << stats.hpdi_high(i) << "\t" << stats.ess(i) << "\t"
                         << stats.rhat(i) << "\n";
        }
    };
    auto write_snp_stats
        = [&snp_names](auto& stream, const auto& stats, Index n)
    {
        for (Index i = 0; i < n; ++i)
        {
            if (i < static_cast<Index>(snp_names.size()))
            {
                stream << snp_names[i];
            }
            else
            {
                stream << "snp" + std::to_string(i + 1);
            }
            stream << "\t" << stats.mean(i) << "\n";
        }
    };

    // saving all parameters except for snp effects
    write_stats(*fixed_, fixed_->size());
    for (const auto& rand : random_)
    {
        write_stats(rand.coeff, rand.coeff.size());
        write_stats(rand.sigma, rand.sigma.size());
    }
    write_stats(residual_, residual_.size());
    write_stats(additive_->variance, additive_->variance.size());
    if (dominant_)
    {
        write_stats(dominant_->variance, dominant_->variance.size());
    }

    // saving snp effects mean
    write_snp_stats(additive_stream, additive_->coeff, additive_->coeff.size());
    if (dominant_)
    {
        write_snp_stats(*dom_stream, dominant_->coeff, dominant_->coeff.size());
    }
}

void MCMCResult::compute_summary_statistics(
    PosteriorSummary& summary,
    const Samples& samples,
    double prob)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return;
    }

    const size_t n_params = samples[0].rows();
    const size_t n_chains = samples.size();
    const size_t n_draws = samples[0].cols();
    int n_threads = Eigen::nbThreads();
    Eigen::setNbThreads(1);

#pragma omp parallel for default(none) \
    shared(samples, prob, summary, n_params, n_chains, n_draws)
    for (size_t r = 0; r < n_params; ++r)
    {
        VectorXd flat_sample(n_chains * n_draws);
        for (size_t chain = 0; chain < n_chains; ++chain)
        {
            flat_sample.segment(chain * n_draws, n_draws)
                = samples[chain].row(r).transpose();
        }

        double mean_val = flat_sample.mean();
        double stddev_val = std::sqrt(
            (flat_sample.array() - mean_val).square().sum()
            / (flat_sample.size() - 1));

        auto [hpdi_low_val, hpdi_high_val] = hpdi(flat_sample, prob);

        summary.mean(r) = mean_val;
        summary.stddev(r) = stddev_val;
        summary.hpdi_low(r) = hpdi_low_val;
        summary.hpdi_high(r) = hpdi_high_val;
    }
    // Restore Eigen threads after all parallel computations
    Eigen::setNbThreads(n_threads);

    summary.ess = effect_sample_size(samples, true);
    summary.rhat = split_gelman_rubin(samples);
}

void MCMCResult::compute_summary_statistics(
    PosteriorSummary& summary,
    const Samples& samples)
{
    if (samples.empty() || samples[0].rows() == 0)
    {
        return;
    }

    const size_t n_params = samples[0].rows();
    const size_t n_chains = samples.size();
    const size_t n_draws = samples[0].cols();
    int n_threads = Eigen::nbThreads();
    Eigen::setNbThreads(1);

#pragma omp parallel for default(none) \
    shared(samples, summary, n_params, n_chains, n_draws)
    for (Index r = 0; r < n_params; ++r)
    {
        VectorXd flat_sample(n_chains * n_draws);
        for (size_t chain = 0; chain < n_chains; ++chain)
        {
            flat_sample.segment(chain * n_draws, n_draws)
                = samples[chain].row(r).transpose();
        }

        double mean_val = flat_sample.mean();
        summary.mean(r) = mean_val;
    }
    Eigen::setNbThreads(n_threads);
}

}  // namespace gelex
