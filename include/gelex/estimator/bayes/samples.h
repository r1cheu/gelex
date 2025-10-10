/**
 * @file samples.h
 * @brief class to save samples in MCMC process easily.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "gelex/estimator/bayes/params.h"

namespace gelex
{

struct MCMCParams;   // Forward declaration
struct BayesStatus;  // Forward declaration
struct BayesModel;   // Forward declaration

// samples for one effect, each element in vector saves one chain's samples. and
// Matrix have the shape (n_params, n_draws)
using Samples = std::vector<Eigen::MatrixXd>;

struct FixedSamples
{
    static FixedSamples create(
        const MCMCParams& params,
        const BayesModel& model);

    Samples coeffs;
    std::vector<std::string> names;
    std::vector<std::string> levels;
    explicit operator bool() const { return !coeffs.empty(); }
};

struct RandomSamples
{
    Samples coeffs;
    Samples effect_variance;
    std::string name;

    void reserve(size_t n_chains)
    {
        coeffs.reserve(n_chains);
        effect_variance.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct AdditiveSamples : RandomSamples
{
    static AdditiveSamples create(
        const MCMCParams& params,
        const BayesModel& model);

    Samples marker_variance;
    void reserve(size_t n_chains)
    {
        RandomSamples::reserve(n_chains);
        marker_variance.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct DominantSamples
{
    static DominantSamples create(
        const MCMCParams& params,
        const BayesModel& model);
    Samples coeffs;
    Samples effect_variance;

    void reserve(size_t n_chains)
    {
        coeffs.reserve(n_chains);
        effect_variance.reserve(n_chains);
    }
    explicit operator bool() const { return !coeffs.empty(); }
};

struct MultiRandomSamples
{
    static MultiRandomSamples create(
        const MCMCParams& params,
        const BayesModel& model);
    std::vector<RandomSamples> chain_samples;
    size_t size() const { return chain_samples.size(); }
    void reserve(size_t n_effects) { chain_samples.reserve(n_effects); }
    explicit operator bool() const { return !chain_samples.empty(); }
};

struct ResidualSamples
{
    static ResidualSamples create(
        const MCMCParams& params,
        const BayesModel& model);
    Samples variance;
    void reserve(size_t n_chains) { variance.reserve(n_chains); }
    explicit operator bool() const { return !variance.empty(); }
};

struct PiSamples
{
    static PiSamples create(const MCMCParams& params, const BayesModel& model);
    Samples prop;
    void reserve(size_t n_chains) { prop.reserve(n_chains); }
    explicit operator bool() const { return !prop.empty(); }
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

   private:
    bool store_fixed_{false};
    bool store_random_{false};
    bool store_pi_{false};
    bool store_dominant_{false};

    FixedSamples fixed_;
    MultiRandomSamples random_;
    AdditiveSamples additive_;
    DominantSamples dominant_;
    ResidualSamples residual_;
    PiSamples pi_;
};
}  // namespace gelex
