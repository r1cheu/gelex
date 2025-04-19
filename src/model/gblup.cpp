#include "gelex/model/gblup.h"

#include <fmt/format.h>

#include <utility>
#include <vector>

#include <armadillo>

namespace gelex
{


GBLUP::GBLUP(std::string formula, dvec phenotype, dmat design_mat_beta)
    : formula_{std::move(formula)},
      phenotype_{std::move(phenotype)},

      design_mat_beta_{std::move(design_mat_beta)}
{
    n_common_effects_ = design_mat_beta_.n_cols;
    beta_ = arma::zeros(n_common_effects_);
}

void GBLUP::add_group_effect(std::string name, sp_dmat design_mat_env)
{
    sp_dmat design_mat{std::move(design_mat_env)};
    group_cov_mats_.emplace_back(design_mat * design_mat.t());
    group_names_.emplace_back(fmt::format("R({})", std::move(name)));
}

void GBLUP::set_model()
{
    n_group_effects_ = group_cov_mats_.size();
    n_individuals_ = phenotype_.n_elem;
    n_genetic_effects_ = genetic_cov_mats_.size();
}

void GBLUP::reset()
{
    group_cov_mats_.clear();
    genetic_cov_mats_.clear();
    genetic_names_.clear();
    group_names_.clear();
}

void GBLUP::add_genetic_effect(std::string name, dmat genetic_covar_mat)
{
    genetic_cov_mats_.emplace_back(std::move(genetic_covar_mat));
    genetic_names_.emplace_back(fmt::format("G({})", std::move(name)));

}

GBLUPParams::GBLUPParams(
    dvec beta,
    dvec sigma,
    dvec proj_y,
    std::vector<std::string> dropped_individuals)
    : beta_{std::move(beta)},
      sigma_{std::move(sigma)},
      proj_y_{std::move(proj_y)},
      dropped_individuals_{std::move(dropped_individuals)} {};
}  // namespace gelex
