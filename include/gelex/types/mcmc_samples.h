/**
 * @file mcmc_samples.h
 * @brief class to save samples in MCMC process easily.
 */

#pragma once

#include <optional>
#include <vector>

#include <Eigen/Core>

// Forward declaration

namespace gelex::bayes
{

struct FixedEffect;
struct RandomEffect;
struct GeneticEffect;
struct AdditiveEffect;
struct DominantEffect;

};  // namespace gelex::bayes

namespace gelex
{

struct MCMCParams;
class BayesState;
class BayesModel;

// samples for one effect, each element in vector saves one chain's samples. and
// Matrix have the shape (n_params, n_draws)
using Samples = std::vector<Eigen::MatrixXd>;
using IntSamples = std::vector<Eigen::MatrixXi>;

struct FixedSamples
{
    FixedSamples(const MCMCParams& params, const bayes::FixedEffect& effect);

    Samples coeffs;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct RandomSamples
{
    RandomSamples(const MCMCParams& params, const bayes::RandomEffect& effect);
    Samples coeffs;
    Samples variance;
    explicit operator bool() const { return !coeffs.empty(); }

   protected:
    RandomSamples(const MCMCParams& params, Eigen::Index n_coeffs);
};

struct BaseMarkerSamples : RandomSamples
{
    BaseMarkerSamples(
        const MCMCParams& params,
        const bayes::GeneticEffect& effect);

    Samples mixture_proportion;
    Samples heritability;
    IntSamples tracker;

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

    Samples variance;
    explicit operator bool() const { return !variance.empty(); }
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
};
}  // namespace gelex
