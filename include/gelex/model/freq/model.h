#pragma once
#include <string>
#include <vector>

#include <armadillo>

#include "gelex/model/freq/effects.h"

namespace gelex
{
/**
 * @class GBLUP
 * @brief Genomic Best Linear Unbiased Prediction (GBLUP) model for genetic
 * analysis.
 *
 * This class handles constructing and managing a GBLUP model, including fixed
 * and random effects, genetic effects, and GxE (genotype by environment).
 *
 * Example usage:
 * ```cpp
 * GBLUP model("y ~ 1 + x", phenotype_data);
 * model.add_fixed_effect({"x"}, {"level1", "level2"}, design_matrix);
 * model.add_genetic_effect("genetic", genetic_matrix, covariance_matrix);
 * ```
 */
class GBLUP
{
   public:
    GBLUP(std::string formula, arma::dvec&& phenotype);

    size_t n_individuals() const { return n_individuals_; }

    size_t n_fixed_effects() const { return fixed_.size(); }
    size_t n_random_effects() const { return random_.n_random_effects(); }
    size_t n_genetic_effects() const { return random_.n_genetic_effects(); }
    size_t n_gxe_effects() const { return random_.n_gxe_effects(); }

    const std::string& formula() const { return formula_; }
    const arma::dvec& phenotype() const { return phenotype_; }

    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        arma::dmat&& design_mat);

    void add_random_effect(std::string&& name, arma::sp_dmat&& design_mat);

    void add_genetic_effect(
        std::string&& name,
        arma::sp_dmat&& design_mat,
        const arma::dmat& cov_mat);

    void add_gxe_effect(
        std::string&& name,
        arma::sp_dmat&& design_mat_genetic,
        const arma::dmat& genetic_cov_mat,
        const arma::dmat& design_mat);

    void add_residual();

    const auto& random() const { return random_; }
    auto& random() { return random_; }

    const auto& fixed() const { return fixed_; }
    auto& fixed() { return fixed_; }

    void clear();

   private:
    std::string formula_;

    size_t n_individuals_{};
    arma::dvec phenotype_;

    RandomEffectManager random_;
    FixedEffect fixed_;
};
}  // namespace gelex
