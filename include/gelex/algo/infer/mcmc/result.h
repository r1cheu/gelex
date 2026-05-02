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

#ifndef GELEX_ALGO_INFER_MCMC_RESULT_H_
#define GELEX_ALGO_INFER_MCMC_RESULT_H_

#include <algorithm>
#include <optional>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/mcmc/samples.h"
#include "gelex/algo/infer/posterior_summary.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

struct AssignmentSummary
{
    explicit AssignmentSummary(const mcmc::AssignmentSamples& samples);
    void compute(const mcmc::AssignmentSamples& sample);

    PosteriorSummary mixture_proportion;
    Eigen::VectorXd pip;
    Eigen::MatrixXd comp_probs;
};

struct ComponentSummary
{
    explicit ComponentSummary(const mcmc::ComponentSamples& samples);
    void compute(const mcmc::ComponentSamples& sample);

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
    explicit GeneticSummary(const mcmc::GeneticSamples& samples);

    void compute(const mcmc::GeneticSamples& sample, double phenotype_var);

    GeneticMode type;
    PosteriorSummary coeffs;
    PosteriorSummary variance;
    PosteriorSummary heritability;
    PosteriorSummary pve;

    std::optional<MixtureSummary> group;
    std::optional<AssignmentSummary> sign;
};

namespace mcmc
{

class Result
{
   public:
    Result(mcmc::Samples&& samples, const BayesModel& model);

    void compute();

    const FixedSummary& fixed() const { return fixed_; }
    const std::vector<RandomSummary>& random() const { return random_; }

    const std::vector<GeneticSummary>& genetics() const { return genetics_; }
    const GeneticSummary* genetic(GeneticMode type) const
    {
        auto it = std::ranges::find(genetics_, type, &GeneticSummary::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const PosteriorSummary& residual() const { return residual_; }

    const Eigen::VectorXd& allele_freq() const { return p_freq_; }

   private:
    mcmc::Samples samples_;

    FixedSummary fixed_;
    std::vector<RandomSummary> random_;
    std::vector<GeneticSummary> genetics_;
    PosteriorSummary residual_;

    double phenotype_var_;

    Eigen::VectorXd p_freq_;
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_RESULT_H_
