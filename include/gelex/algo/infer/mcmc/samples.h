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

#ifndef GELEX_ALGO_INFER_MCMC_SAMPLES_H_
#define GELEX_ALGO_INFER_MCMC_SAMPLES_H_

#include <algorithm>
#include <cstddef>
#include <memory>
#include <optional>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "gelex/algo/infer/fixed_samples.h"
#include "gelex/bayes/genetic/prior_state.h"
#include "gelex/bayes/model.h"
#include "gelex/bayes/prior.h"
#include "gelex/bayes/state.h"
#include "gelex/infra/stats/detail/running_stats.h"
#include "gelex/infra/stats/result.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

struct GeneticDesign;
struct ResidualState;

};  // namespace gelex::bayes

namespace gelex
{

class BayesModel;

namespace mcmc
{

class Writer;

struct AssignmentSamples
{
    AssignmentSamples(
        Eigen::Index n_snps,
        Eigen::Index n_proportions,
        bool estimate_pi);

    void store(const Eigen::VectorXi& assignment);
    void store(const bayes::MixtureState& state);

    auto n_snps() const -> Eigen::Index { return comp_counts_.rows(); }
    auto n_proportions() const -> Eigen::Index { return comp_counts_.cols(); }
    auto estimate_pi() const -> bool { return estimate_pi_; }
    auto proportion() const -> stats::RunningStatsResult
    {
        return proportion_stats_.result();
    }
    auto component_probs() const -> Eigen::MatrixXd
    {
        return comp_counts_ / static_cast<double>(n_samples_);
    }

   private:
    bool estimate_pi_;
    stats::detail::RunningStats proportion_stats_;
    Eigen::MatrixXd comp_counts_;
    std::size_t n_samples_{0};
};

struct ComponentSamples
{
    ComponentSamples(
        Eigen::Index n_snps,
        Eigen::Index n_proportions,
        bool estimate_pi);

    void store(
        const bayes::ComponentState& component,
        const Eigen::VectorXi& mixture_assignment);
    void store(
        const bayes::ComponentState& component,
        const bayes::MixtureState& mixture);

    auto comp_var() const -> stats::RunningStatsResult
    {
        return comp_var_stats_.result();
    }

    AssignmentSamples assignment;

   private:
    stats::detail::RunningStats comp_var_stats_;
};

using MixtureSamples = std::variant<AssignmentSamples, ComponentSamples>;

inline auto assignment(const MixtureSamples& s) -> const AssignmentSamples&
{
    return std::visit(
        []<typename T>(const T& v) -> const AssignmentSamples&
        {
            if constexpr (std::is_same_v<T, AssignmentSamples>)
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

struct GeneticSamples
{
    GeneticSamples(const bayes::GeneticDesign& design, GeneticMode mode);
    void store(const BayesState& state);

    auto n_coeffs() const -> Eigen::Index { return n_coeffs_; }
    auto coeffs() const -> stats::RunningStatsResult
    {
        return coeffs_stats_.result();
    }
    auto variance() const -> stats::RunningStatsResult
    {
        return variance_stats_.result();
    }
    auto heritability() const -> stats::RunningStatsResult
    {
        return heritability_stats_.result();
    }

    GeneticMode type;
    std::optional<MixtureSamples> group;
    std::optional<AssignmentSamples> sign;

   private:
    Eigen::Index n_coeffs_;
    stats::detail::RunningStats coeffs_stats_;
    stats::detail::RunningStats variance_stats_;
    stats::detail::RunningStats heritability_stats_;
};

struct ResidualSamples
{
    ResidualSamples() = default;
    void store(const bayes::ResidualState& state);

    auto variance() const -> stats::RunningStatsResult
    {
        return variance_stats_.result();
    }

   private:
    stats::detail::RunningStats variance_stats_;
};

class Samples
{
   public:
    Samples(const Samples&) = delete;
    auto operator=(const Samples&) -> Samples& = delete;
    Samples(Samples&&) noexcept;
    auto operator=(Samples&&) noexcept -> Samples&;
    ~Samples();

    Samples(
        const BayesModel& model,
        const bayes::BayesPrior& prior,
        const BayesState& state,
        std::string_view sample_prefix,
        Eigen::Index n_records);
    void store(const BayesState& states);
    void finalize();

    const FixedSamples& fixed() const { return fixed_; }
    const std::vector<RandomSamples>& random() const { return random_; }

    const std::vector<GeneticSamples>& genetics() const { return genetics_; }
    const GeneticSamples* genetic(GeneticMode type) const
    {
        auto it = std::ranges::find(genetics_, type, &GeneticSamples::type);
        return it != genetics_.end() ? &*it : nullptr;
    }

    const ResidualSamples& residual() const { return residual_; }

   private:
    FixedSamples fixed_;
    std::vector<RandomSamples> random_;
    std::vector<GeneticSamples> genetics_;
    ResidualSamples residual_;
    std::unique_ptr<mcmc::Writer> writer_;
};

}  // namespace mcmc

}  // namespace gelex

#endif  // GELEX_ALGO_INFER_MCMC_SAMPLES_H_
