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
        arma::dmat&& design_matrix);
    void add_random_effect(std::string&& name, arma::dmat&& design_matrix);
    void add_genetic_effect(
        std::string&& name,
        arma::dmat&& genotype,
        BayesAlphabet type);

    const std::string& formula() const { return formula_; }

    const auto& fixed() const { return *fixed_; }
    const auto& random() const { return random_; }
    const auto& genetic() const { return genetic_; }
    const auto& residual() const { return residual_; }

    auto& fixed() { return *fixed_; }
    auto& random() { return random_; }
    auto& genetic() { return genetic_; }
    auto& residual() { return residual_; }

    // Existing method (keep for expert users)
    void set_sigma_prior_manul(
        const std::string& name,
        double nu,
        double s2,
        double init_sigma);

    // New overload with default calculation
    void set_sigma_prior(
        const std::string& name,
        double sigma_prop,
        double nu = 4.0);

    // Keep existing pi prior method
    void set_pi_prior(const std::string& name, const arma::dvec& pi);

    std::string prior_summary() const;

    const arma::dvec& phenotype() const { return phenotype_; }
    double phenotype_var() const { return phenotype_var_; }
    size_t n_individuals() const { return n_individuals_; }

   private:
    std::string formula_;

    size_t n_individuals_{};

    arma::dvec phenotype_;
    double phenotype_var_{};

    std::unique_ptr<bayes::FixedEffect> fixed_;
    bayes::RandomEffectManager random_;
    bayes::GeneticEffectManager genetic_;
    bayes::Residual residual_;
};

struct BayesStatus
{
    explicit BayesStatus(const BayesModel& model)
        : fixed(bayes::FixedEffectState(model.fixed().design_matrix.n_cols)),
          random(create_thread_states(model.random())),
          genetic(create_thread_states(model.genetic())),
          residual(model.residual())
    {
    }

    void compute_heritability();
    bayes::FixedEffectState fixed;
    std::vector<bayes::RandomEffectState> random;
    std::vector<bayes::GeneticEffectState> genetic;
    bayes::Residual residual;
};

}  // namespace gelex
