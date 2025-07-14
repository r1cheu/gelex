#pragma once
#include <armadillo>
#include <vector>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/model.h"

namespace gelex
{

struct MCMCParams;  // Forward declaration

struct RandomGroup
{
    std::vector<arma::dcube> coeffs;
    std::vector<arma::dcube> sigmas;
    size_t size() const { return coeffs.size(); }
};

struct GeneticGroup : RandomGroup
{
    std::vector<arma::dcube> genetic_var;
    std::vector<arma::dcube> heritability;
    std::vector<arma::dcube> pi;
    size_t size() const { return coeffs.size(); }
};

class MCMCSamples
{
   public:
    explicit MCMCSamples(const MCMCParams& params, const BayesModel& model);
    void store(const BayesStatus& status, size_t record_idx, size_t chain_idx);

    const auto& fixed() const { return fixed_; }
    const auto& random() const { return random_; }
    const auto& genetic() const { return genetic_; }
    const auto& residual() const { return residual_; }

   private:
    size_t n_records_;
    size_t n_chains_;

    arma::dcube fixed_;
    RandomGroup random_;
    GeneticGroup genetic_;
    arma::dcube residual_;

    void init_group(
        RandomGroup& group,
        const bayes::RandomEffectManager& effects) const;

    void init_group(
        GeneticGroup& group,
        const bayes::GeneticEffectManager& effects) const;

    void store_group(
        const std::vector<bayes::RandomEffectState>& status,
        size_t record_idx,
        size_t chain_idx);

    void store_group(
        const std::vector<bayes::GeneticEffectState>& status,
        size_t record_idx,
        size_t chain_idx);
};
}  // namespace gelex
