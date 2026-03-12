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

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/types/genetic_effect_type.h"
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
    explicit PosteriorSummary(RunningStatsResult result);
    PosteriorSummary() = default;

    Eigen::Index size() const { return mean.size(); }

    Eigen::VectorXd mean;
    Eigen::VectorXd stddev;
};

struct FixedSummary
{
    explicit FixedSummary(const FixedSamples& sample)
        : names(sample.names), levels(sample.levels), coeffs(sample.n_coeffs())
    {
    }

    void compute(const FixedSamples& sample);

    template <typename F>
    void for_each_term(const F& fn) const
    {
        Eigen::Index coeff_idx = 0;
        for (size_t group_idx = 0; group_idx < levels.size(); ++group_idx)
        {
            const auto& group_levels = levels[group_idx];
            if (group_levels)
            {
                for (const auto& level : *group_levels)
                {
                    fn(names[group_idx] + "_" + level, coeff_idx);
                    ++coeff_idx;
                }
            }
            else
            {
                fn(names[group_idx], coeff_idx);
                ++coeff_idx;
            }
        }
    }

    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    PosteriorSummary coeffs;
};

struct RandomSummary
{
    explicit RandomSummary(const RandomSamples& sample)
        : name(sample.name),
          levels(sample.levels),
          coeffs(sample.n_coeffs()),
          variance(1)
    {
    }

    void compute(const RandomSamples& sample);

    template <typename F>
    void for_each_term(const F& fn) const
    {
        Eigen::Index coeff_idx = 0;
        if (levels)
        {
            for (const auto& level : *levels)
            {
                fn(name.empty() ? level : name + "_" + level, coeff_idx);
                ++coeff_idx;
            }
        }
        else
        {
            fn(name, coeff_idx);
            ++coeff_idx;
        }
    }

    std::string name;
    std::optional<std::vector<std::string>> levels;
    PosteriorSummary coeffs;
    PosteriorSummary variance;
};

struct MixtureSummary
{
    explicit MixtureSummary(const MixtureSamples& samples)
        : pip(Eigen::VectorXd::Zero(samples.n_snps())),
          comp_probs(
              Eigen::MatrixXd::Zero(samples.n_snps(), samples.n_proportions()))
    {
        if (samples.estimate_pi())
        {
            mixture_proportion = PosteriorSummary(samples.n_proportions());
        }

        if (samples.n_proportions() > 2)
        {
            component_variance = PosteriorSummary(samples.n_proportions() - 1);
        }
    }

    void compute(const MixtureSamples& sample);

    PosteriorSummary mixture_proportion;
    PosteriorSummary component_variance;
    Eigen::VectorXd pip;
    Eigen::MatrixXd comp_probs;
};

struct SignSummary
{
    explicit SignSummary(const SignSamples& samples);
    void compute(const SignSamples& samples);

    PosteriorSummary positive_prob;
};

struct GeneticSummary
{
    explicit GeneticSummary(const GeneticSamples& samples)
        : type(samples.type),
          coeffs(samples.n_coeffs()),
          variance(1),
          heritability(1),
          pve(samples.n_coeffs())
    {
        if (samples.mixture)
        {
            mixture.emplace(*samples.mixture);
        }
        if (samples.sign)
        {
            sign.emplace(*samples.sign);
        }
    }

    void compute(const GeneticSamples& sample, double phenotype_var);

    GeneticEffectType type;
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary heritability;
    PosteriorSummary pve;

    std::optional<MixtureSummary> mixture;
    std::optional<SignSummary> sign;
};

class MCMCResult
{
   public:
    MCMCResult(MCMCSamples&& samples, const BayesModel& model);

    void compute();

    const FixedSummary& fixed() const { return fixed_; }
    const std::vector<RandomSummary>& random() const { return random_; }

    const std::vector<GeneticSummary>& genetics() const { return genetics_; }
    const GeneticSummary* genetic(GeneticEffectType type) const
    {
        auto it = std::ranges::find(genetics_, type, &GeneticSummary::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const PosteriorSummary& residual() const { return residual_; }

    const Eigen::VectorXd& allele_freq() const { return p_freq_; }

   private:
    MCMCSamples samples_;

    FixedSummary fixed_;
    std::vector<RandomSummary> random_;
    std::vector<GeneticSummary> genetics_;
    PosteriorSummary residual_;

    double phenotype_var_;

    Eigen::VectorXd p_freq_;
};

}  // namespace gelex

#endif  // GELEX_TYPES_MCMC_RESULTS_H_
