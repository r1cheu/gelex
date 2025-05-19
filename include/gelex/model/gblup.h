#pragma once
#include <string>
#include <vector>
#include "gelex/model/effects.h"

#include <armadillo>

namespace gelex
{
using arma::dcube;
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
using arma::uvec;
using arma::uword;
using std::string;

class GBLUP
{
   public:
    GBLUP(string formula, dvec phenotype);

    size_t n_individuals() const { return n_individuals_; }

    size_t n_fixed_effects() const { return fixed_.size(); }
    size_t n_random_effects() const { return random_.n_random_effects(); }
    size_t n_genetic_effects() const { return random_.n_genetic_effects(); }
    size_t n_gxe_effects() const { return random_.n_gxe_effects(); }

    const string& formula() const { return formula_; }
    const dvec& phenotype() const { return phenotype_; }

    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        dmat&& design_mat_beta);

    void add_random_effect(string&& name, sp_dmat&& design_mat_random);

    void add_genetic_effect(
        string&& name,
        sp_dmat&& design_mat_genetic,
        const dmat& genetic_covar_mat);

    void add_gxe_effect(
        string&& name,
        sp_dmat&& design_mat_genetic,
        const dmat& genetic_cov_mat,
        const dmat& design_mat_group);

    void add_residual();

    const RandomEffectManager& random() const { return random_; }
    RandomEffectManager& random() { return random_; }

    const FixedEffect& fixed() const { return fixed_; }
    FixedEffect& fixed() { return fixed_; }

    void clear();

   private:
    std::string formula_;

    size_t n_individuals_{};
    dvec phenotype_;

    RandomEffectManager random_;
    FixedEffect fixed_;
};

struct GBLUPParams
{
    dvec beta;
    dvec sigma;
    dvec proj_y;
    std::vector<std::string> dropped_ids;
};

}  // namespace gelex
