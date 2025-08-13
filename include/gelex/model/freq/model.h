#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include <armadillo>

#include "gelex/model/freq/effects.h"

namespace gelex
{

/**
 * @class GBLUPEffect
 * @brief add struct for representing an effect in the GBLUP model. Not used by
 * user.
 */
struct GBLUPEffect
{
    std::string name;
    effect_type type;
    double sigma;
    double se;
};

using TotalEffects = Effects<GBLUPEffect>;

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
    size_t n_fixed_effects() const { return fixed_->size(); }
    size_t n_random_effects() const { return random_.size(); }
    size_t n_genetic_effects() const { return genetic_.size(); }
    size_t n_gxe_effects() const { return gxe_.size(); }

    const std::string& formula() const { return formula_; }
    const arma::dvec& phenotype() const { return phenotype_; }
    arma::dvec var_comp() const;
    size_t n_var_comp() const;

    std::unordered_map<std::string, arma::dvec> u()
    {
        std::unordered_map<std::string, arma::dvec> u;
        for (const auto& effect : random_)
        {
            u[effect.name] = effect.coeff;
        }

        for (const auto& effect : genetic_)
        {
            u[effect.name] = effect.coeff;
        }
        return u;
    }

    const FixedEffect& fixed() const { return *fixed_; }
    const Effects<GBLUPEffect>& effects() const { return effects_; }

    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        arma::dmat&& design_matrix);
    void add_random_effect(std::string&& name, arma::sp_dmat&& design_matrix);
    void add_genetic_effect(
        std::string&& name,
        arma::sp_dmat&& design_matrix,
        const arma::dmat& genetic_relationship_matrix);
    void add_gxe_effect(
        std::string&& name,
        arma::sp_dmat&& genetic_design_matrix,
        const arma::dmat& genetic_relationship_matrix,
        arma::sp_dmat&& env_design_matrix);
    void clear();

    friend class Optimizer;
    friend class Estimator;
    friend class GBLUPPredictor;

   private:
    void set_var_comp(const arma::dvec& var_comp);
    void init_var_comp(double var);

    std::string formula_;

    size_t n_individuals_{};
    arma::dvec phenotype_;
    Effects<GBLUPEffect> effects_;

    std::unique_ptr<FixedEffect> fixed_;
    RandomEffects random_;
    GeneticEffects genetic_;
    GxEEffects gxe_;
    RandomEffect residual_{"e", arma::sp_dmat()};
};

}  // namespace gelex
