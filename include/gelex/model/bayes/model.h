#pragma once

#include <string>
#include <vector>

#include <armadillo>

#include "gelex/model/bayes/effects.h"

namespace gelex
{
/**
 * @class BayesModel
 * @brief Bayesian model for genetic analysis with fixed, random, and genetic
 * effects.
 *
 * Supports adding and managing effects (fixed, random, genetic) and provides
 * accessors for model components (mu, residuals, etc.). Initialized with a
 * formula and phenotype.
 */
class BayesModel
{
   public:
    /**
     * @brief
     *
     * @param formula formula string, only for showing the model structure.
     * @param phenotype moveable arma dvec.
     */
    BayesModel(std::string formula, arma::dvec&& phenotype);
    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        arma::dmat&& design_mat);
    void add_random_effect(std::string&& name, arma::dmat&& design_mat);
    void add_genetic_effect(
        std::string&& name,
        arma::dmat&& genotype,
        BayesAlphabet type);

    const std::string& formula() const { return formula_; }

    const auto& mu() const { return mu_; }
    const auto& fixed() const { return fixed_; }
    const auto& random() const { return random_; }
    const auto& genetic() const { return genetic_; }
    const auto& residual() const { return residual_; }

    auto& mu() { return mu_; }
    auto& fixed() { return fixed_; }
    auto& random() { return random_; }
    auto& genetic() { return genetic_; }
    auto& residual() { return residual_; }

    void set_sigma_prior(const std::string& name, double nu, double s2);
    void set_pi_prior(const std::string& name, const arma::dvec& pi);
    void set_residual_prior(double nu, double s2);

    const arma::dvec& phenotype() const { return phenotype_; }
    double phenotype_var() const { return phenotype_var_; }

    size_t n_individuals() const { return n_individuals_; }

   private:
    std::string formula_;

    size_t n_individuals_{};

    arma::dvec phenotype_;
    double phenotype_var_{};

    Mu mu_;
    std::unique_ptr<FixedEffectDesign> fixed_;
    RandomEffectDesignManager random_;
    GeneticEffectDesignManager genetic_;
    Residual residual_;
};

struct BayesStatus
{
    explicit BayesStatus(const BayesModel& model)
        : mu(model.mu()),
          fixed(
              model.fixed() ? FixedEffectState(model.fixed()->design_mat.n_cols)
                            : FixedEffectState(0)),
          random(create_thread_states(model.random())),
          genetic(create_thread_states(model.genetic())),
          residual(model.residual())
    {
    }
    Mu mu;
    FixedEffectState fixed;
    std::vector<RandomEffectState> random;
    std::vector<GeneticEffectState> genetic;
    Residual residual;
};

}  // namespace gelex
