/**
 * @file samples.h
 * @brief class to save samples in MCMC process easily.
 */

#pragma once

#include <optional>
#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/estimator/bayes/params.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration
class BayesState;   // Forward declaration
class BayesModel;   // Forward declaration

// samples for one effect, each element in vector saves one chain's samples. and
// Matrix have the shape (n_params, n_draws)
using Samples = std::vector<Eigen::MatrixXd>;

struct FixedSamples
{
    FixedSamples(const MCMCParams& params, const bayes::FixedEffect& effect);

    Samples coeffs;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct RandomSamples
{
    RandomSamples(const MCMCParams& params, const bayes::RandomEffect& effect)
        : RandomSamples(params, effect.design_matrix.cols())
    {
    }
    Samples coeffs;
    Samples variance;
    explicit operator bool() const { return !coeffs.empty(); }

   protected:
    RandomSamples(const MCMCParams& params, Eigen::Index n_coeffs);
};

struct AdditiveSamples : RandomSamples
{
    AdditiveSamples(
        const MCMCParams& params,
        const bayes::AdditiveEffect& state)
        : RandomSamples(params, bayes::get_cols(state.design_matrix))
    {
    }
};

struct DominantSamples : RandomSamples
{
    DominantSamples(
        const MCMCParams& params,
        const bayes::DominantEffect& effect)
        : RandomSamples(params, bayes::get_cols(effect.design_matrix))
    {
        ratios.reserve(params.n_chains);
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            ratios.emplace_back(
                bayes::get_cols(effect.design_matrix), params.n_records);
        }
    }
    Samples ratios;
};

struct ResidualSamples
{
    explicit ResidualSamples(const MCMCParams& params)
    {
        variance.reserve(params.n_chains);
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            variance.emplace_back(1, params.n_records);
        }
    }
    Samples variance;
    explicit operator bool() const { return !variance.empty(); }
};

struct PiSamples
{
    PiSamples(const MCMCParams& params, const bayes::AdditiveEffect& effect)
    {
        prop.reserve(params.n_chains);
        const Eigen::Index n_props = effect.pi.size();
        for (Eigen::Index i = 0; i < params.n_chains; ++i)
        {
            prop.emplace_back(n_props, params.n_records);
        }
    }

    Samples prop;
    explicit operator bool() const { return !prop.empty(); }
};

class MCMCSamples
{
   public:
    MCMCSamples(const MCMCParams& params, const BayesModel& model);
    void store(
        const BayesState& states,
        Eigen::Index record_idx,
        Eigen::Index chain_idx);

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
    std::optional<PiSamples> pi_;
};
}  // namespace gelex
