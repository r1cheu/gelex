#pragma once

#include <string>
#include <vector>

#include <armadillo>

#include "gelex/model/effects/bayes_effects.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;

class Bayes
{
   public:
    Bayes(std::string formula, dvec&& phenotype);
    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        dmat&& design_mat);
    void add_random_effect(std::string&& name, dmat&& design_mat);
    void
    add_genetic_effect(std::string&& name, dmat&& genotype, BayesAlphabet type);

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

    const dvec& phenotype() const { return phenotype_; }
    double phenotype_var() const { return phenotype_var_; }

    size_t n_individuals() const { return n_individuals_; }

   private:
    std::string formula_;

    size_t n_individuals_{};

    dvec phenotype_;
    double phenotype_var_{};

    bayes::Mu mu_;
    bayes::FixedEffect fixed_;
    bayes::RandomEffectManager random_;
    bayes::GeneticEffectManager genetic_;
    bayes::Residual residual_;

    arma::rowvec allele_freq_;
};

}  // namespace gelex
