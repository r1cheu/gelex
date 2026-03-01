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

/**
 * @file mcmc_samples.h
 * @brief class to save samples in MCMC process easily.
 */

#ifndef GELEX_TYPES_MCMC_SAMPLES_H_
#define GELEX_TYPES_MCMC_SAMPLES_H_

#include <memory>
#include <optional>
#include <vector>

#include <Eigen/Core>

// Forward declaration

namespace gelex::detail
{

template <typename eT>
class BinaryWriter;

}  // namespace gelex::detail

namespace gelex::bayes
{

struct RandomEffect;
struct GeneticEffect;
struct AdditiveEffect;
struct DominantEffect;

};  // namespace gelex::bayes

namespace gelex
{

struct FixedEffect;
struct MCMCParams;
class BayesState;
class BayesModel;

struct FixedSamples
{
    FixedSamples(const MCMCParams& params, const FixedEffect& effect);

    Eigen::MatrixXd coeffs;
    explicit operator bool() const { return coeffs.size() > 0; }
};

struct RandomSamples
{
    RandomSamples(const MCMCParams& params, const bayes::RandomEffect& effect);
    Eigen::MatrixXd coeffs;
    Eigen::RowVectorXd variance;
    explicit operator bool() const { return coeffs.size() > 0; }

   protected:
    RandomSamples(const MCMCParams& params, Eigen::Index n_coeffs);
};

struct BaseMarkerSamples : RandomSamples
{
    BaseMarkerSamples(
        const MCMCParams& params,
        const bayes::GeneticEffect& effect);

    Eigen::MatrixXd mixture_proportion;
    Eigen::RowVectorXd heritability;
    Eigen::MatrixXi tracker;
    Eigen::MatrixXd component_variance;

    Eigen::Index n_proportions
        = 0;  // load the number of prop for no-estimate-pi models.
};

struct AdditiveSamples : BaseMarkerSamples
{
    AdditiveSamples(
        const MCMCParams& params,
        const bayes::AdditiveEffect& effect);
};

struct DominantSamples : BaseMarkerSamples
{
    DominantSamples(
        const MCMCParams& params,
        const bayes::DominantEffect& effect);
};

struct ResidualSamples
{
    explicit ResidualSamples(const MCMCParams& params);

    Eigen::RowVectorXd variance;
    explicit operator bool() const { return variance.size() > 0; }
};

class MCMCSamples
{
   public:
    MCMCSamples(const MCMCSamples&) = delete;
    auto operator=(const MCMCSamples&) -> MCMCSamples& = delete;
    MCMCSamples(MCMCSamples&&) noexcept;
    auto operator=(MCMCSamples&&) noexcept -> MCMCSamples&;
    ~MCMCSamples();

    MCMCSamples(
        const MCMCParams& params,
        const BayesModel& model,
        std::string_view sample_prefix);
    void store(const BayesState& states, Eigen::Index record_idx);

    const FixedSamples* fixed() const
    {
        return fixed_ ? &fixed_.value() : nullptr;
    }
    const std::vector<RandomSamples>& random() const { return random_; }
    const AdditiveSamples* additive() const
    {
        return additive_ ? &additive_.value() : nullptr;
    }
    const DominantSamples* dominant() const
    {
        return dominant_ ? &dominant_.value() : nullptr;
    }
    const ResidualSamples& residual() const { return residual_; }

   private:
    std::optional<FixedSamples> fixed_;
    std::vector<RandomSamples> random_;
    std::optional<AdditiveSamples> additive_;
    std::optional<DominantSamples> dominant_;
    ResidualSamples residual_;
    std::unique_ptr<detail::BinaryWriter<double>> add_writer_;
    std::unique_ptr<detail::BinaryWriter<double>> dom_writer_;
    std::unique_ptr<detail::BinaryWriter<double>> scalar_writer_;
};
}  // namespace gelex

#endif  // GELEX_TYPES_MCMC_SAMPLES_H_
