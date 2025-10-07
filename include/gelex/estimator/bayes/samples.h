/**
 * @file samples.h
 * @brief class to save samples in MCMC process easily.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration

// samples for one effect, each element in vector saves one chain's samples. and
// Matrix have the shape (n_params, n_draws)
using Samples = std::vector<Eigen::MatrixXd>;

struct FixedSamples
{
    Samples coeff;
    std::vector<std::string> names;
    std::vector<std::string> levels;
    explicit operator bool() const { return !coeff.empty(); }
};

struct RandomSamples
{
    Samples coeffs;
    Samples sigmas;
    std::string name;

    void reserve(size_t n_chains)
    {
        coeffs.reserve(n_chains);
        sigmas.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct AdditiveSamples : RandomSamples
{
    Samples variance;
    void reserve(size_t n_chains)
    {
        RandomSamples::reserve(n_chains);
        variance.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct DominantSamples
{
    Samples coeffs;
    Samples variance;
    void reserve(size_t n_chains)
    {
        coeffs.reserve(n_chains);
        variance.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct MultiRandomSamples
{
    std::vector<RandomSamples> chain_samples;
    size_t size() const { return chain_samples.size(); }
    void reserve(size_t n_effects) { chain_samples.reserve(n_effects); }
    explicit operator bool() const { return !chain_samples.empty(); }
};

class MCMCSamples
{
   public:
    explicit MCMCSamples(const MCMCParams& params, const BayesModel& model);
    void store(
        const BayesStatus& status,
        Eigen::Index record_idx,
        Eigen::Index chain_idx);

    const auto& fixed() const { return fixed_; }
    const auto& random() const { return random_; }
    const auto& additive() const { return additive_; }
    const auto& dominant() const { return dominant_; }
    const auto& residual() const { return residual_; }
    const auto& bim_file_path() const { return bim_file_path_; }

   private:
    size_t n_records_;
    size_t n_chains_;

    bool store_pi_{false};
    std::string bim_file_path_;

    FixedSamples fixed_;
    MultiRandomSamples random_;
    AdditiveSamples additive_;
    DominantSamples dominant_;
    Samples residual_;
    Samples pi_;

    void init_random(const bayes::RandomEffectManager& effects);

    void store_random(
        const std::vector<bayes::RandomStatus>& status,
        Eigen::Index record_idx,
        Eigen::Index chain_idx);
};
}  // namespace gelex
