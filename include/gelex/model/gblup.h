#pragma once
#include <sys/types.h>
#include <cstdint>
#include <string>
#include <utility>
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

    uint64_t n_individuals() const { return n_individuals_; }
    uint64_t n_common_effects() const { return n_common_effects_; }
    uint64_t n_group_effects() const { return effects_.n_group_effects(); }
    uint64_t n_genetic_effects() const { return effects_.n_genetic_effects(); }
    uint64_t n_gxe_effects() const { return effects_.n_gxe_effects(); }

    const string& formula() const { return formula_; }
    const dvec& phenotype() const { return phenotype_; }
    const dmat& design_mat_beta() const { return design_mat_beta_; }
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }

    const std::vector<std::string>& fixed_effect_names() const
    {
        return fixed_eff_names_;
    }

    const std::vector<std::string>& fixed_effect_levles() const
    {
        return fixed_eff_levels_;
    }

    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        dmat&& design_mat_beta);
    void add_group_effect(string name, sp_dmat design_mat_group);
    void add_genetic_effect(
        string name,
        sp_dmat design_mat_genetic,
        const dmat& genetic_covar_mat);
    void add_gxe_effect(
        string name,
        sp_dmat design_mat_genetic,
        const dmat& genetic_cov_mat,
        const dmat& design_mat_group);

    const EffectManager& effect() const { return effects_; }
    EffectManager& effect() { return effects_; }

    void set_sigma(dvec sigma)
    {
        sigma_ = std::move(sigma);
        uint64_t idx{};
        for (auto& eff : effects_)
        {
            eff.sigma = sigma_.at(idx);
            ++idx;
        }
        effects_.set_residual(sigma_.back());
    }

    void set_beta(dvec beta) { beta_ = std::move(beta); }

    void clear();

   private:
    uint64_t n_individuals_{};
    uint64_t n_common_effects_{};

    std::string formula_;
    std::vector<std::string> fixed_eff_names_;
    std::vector<std::string> fixed_eff_levels_;
    dvec phenotype_;

    dmat design_mat_beta_;

    EffectManager effects_;

    dvec beta_;
    dvec sigma_;
};

struct GBLUPParams
{
    dvec beta;
    dvec sigma;
    dvec proj_y;
    std::vector<std::string> dropped_ids;
};

}  // namespace gelex
