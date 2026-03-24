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
#include <variant>
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

struct AssignmentSummary
{
    explicit AssignmentSummary(const AssignmentSamples& samples);
    void compute(const AssignmentSamples& sample);

    PosteriorSummary mixture_proportion;
    Eigen::VectorXd pip;
    Eigen::MatrixXd comp_probs;
};

struct ComponentSummary
{
    explicit ComponentSummary(const ComponentSamples& samples);
    void compute(const ComponentSamples& sample);

    AssignmentSummary assignment;
    PosteriorSummary component_variance;
};

using MixtureSummary = std::variant<AssignmentSummary, ComponentSummary>;

inline auto assignment(const MixtureSummary& s) -> const AssignmentSummary&
{
    return std::visit(
        []<typename T>(const T& v) -> const AssignmentSummary&
        {
            if constexpr (std::is_same_v<T, AssignmentSummary>)
            {
                return v;
            }
            else
            {
                return v.assignment;
            }
        },
        s);
}

struct GeneticSummary
{
    explicit GeneticSummary(const GeneticSamples& samples);

    void compute(const GeneticSamples& sample, double phenotype_var);

    GeneticKind type;
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary heritability;
    PosteriorSummary pve;

    std::optional<MixtureSummary> group;
    std::optional<AssignmentSummary> sign;
};

class MCMCResult
{
   public:
    MCMCResult(MCMCSamples&& samples, const BayesModel& model);

    void compute();

    const FixedSummary& fixed() const { return fixed_; }
    const std::vector<RandomSummary>& random() const { return random_; }

    const std::vector<GeneticSummary>& genetics() const { return genetics_; }
    const GeneticSummary* genetic(GeneticKind type) const
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
