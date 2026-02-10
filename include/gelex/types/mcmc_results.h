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

#ifndef GELEX_TYPES_MCMC_RESULTS_H_
#define GELEX_TYPES_MCMC_RESULTS_H_

#include <optional>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/mcmc_samples.h"

namespace gelex
{

struct PosteriorSummary
{
    explicit PosteriorSummary(Eigen::Index n_params)
        : mean(Eigen::VectorXd::Zero(n_params)),
          stddev(Eigen::VectorXd::Zero(n_params))
    {
    }
    PosteriorSummary() = default;

    Eigen::Index size() const { return mean.size(); }

    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
};

struct FixedSummary
{
    explicit FixedSummary(const FixedSamples& sample)
        : coeffs(sample.coeffs[0].rows())
    {
    }

    PosteriorSummary coeffs;
};

struct RandomSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : coeffs(sample.coeffs[0].rows()), variance(1)
    {
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

struct BaseMarkerSummary
{
    explicit BaseMarkerSummary(const BaseMarkerSamples& samples)
        : coeffs(samples.coeffs[0].rows()),
          variance(1),
          heritability(1),
          pve(samples.coeffs[0].rows())
    {
        if (!samples.tracker.empty())  // mixture model
        {
            pip = Eigen::VectorXd::Zero(samples.tracker[0].rows());
            comp_probs = Eigen::MatrixXd::Zero(
                samples.tracker[0].rows(), samples.n_proportions);
        }

        if (!samples.mixture_proportion.empty())
        {
            mixture_proportion
                = PosteriorSummary(samples.mixture_proportion[0].rows());
        }

        if (!samples.component_variance.empty())
        {
            component_variance
                = PosteriorSummary(samples.component_variance[0].rows());
        }
    }

    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary heritability;
    PosteriorSummary pve;

    PosteriorSummary mixture_proportion;
    PosteriorSummary component_variance;
    Eigen::VectorXd pip;         // Posterior inclusion probability
    Eigen::MatrixXd comp_probs;  // Per-component posterior probabilities
};

struct AdditiveSummary : BaseMarkerSummary
{
    explicit AdditiveSummary(const AdditiveSamples& samples)
        : BaseMarkerSummary(samples)
    {
    }
};

struct DominantSummary : BaseMarkerSummary
{
    explicit DominantSummary(const DominantSamples& samples)
        : BaseMarkerSummary(samples)
    {
    }
};

class MCMCResult
{
   public:
    explicit MCMCResult(
        MCMCSamples&& samples,
        const BayesModel& model,
        double prob = 0.9);

    /**
     * @brief Compute posterior statistics.
     *
     * If prob is provided, uses it as the probability threshold for
     * computation. Otherwise, uses default prob.
     *
     * @param prob Optional probability threshold for computation.
     */
    void compute(std::optional<double> prob = std::nullopt);

    const FixedSummary* fixed() const
    {
        return fixed_ ? &fixed_.value() : nullptr;
    }
    const std::vector<RandomSummary>& random() const { return random_; }
    const AdditiveSummary* additive() const
    {
        return additive_ ? &additive_.value() : nullptr;
    }
    const DominantSummary* dominant() const
    {
        return dominant_ ? &dominant_.value() : nullptr;
    }
    const PosteriorSummary& residual() const { return residual_; }

   private:
    friend class SnpEffectsWriter;
    friend class SnpQuantGeneticWriter;

    MCMCSamples samples_;

    std::optional<FixedSummary> fixed_;
    std::vector<RandomSummary> random_;
    std::optional<AdditiveSummary> additive_;
    std::optional<DominantSummary> dominant_;
    PosteriorSummary residual_;

    double prob_;
    double phenotype_var_;

    Eigen::VectorXd p_freq;
};

}  // namespace gelex

#endif  // GELEX_TYPES_MCMC_RESULTS_H_
