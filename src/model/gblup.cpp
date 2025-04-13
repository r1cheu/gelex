#include "gelex/model/gblup.h"

#include <utility>
#include <vector>

#include <armadillo>

namespace gelex
{

GBLUP::GBLUP(dvec phenotype, dmat design_mat_beta)
    : phenotype_{std::move(phenotype)},
      design_mat_beta_{std::move(design_mat_beta)}
{
    n_common_effects_ = design_mat_beta_.n_cols;
    beta_ = arma::zeros(n_common_effects_);
}

void GBLUP::add_group_effect(std::string name, sp_dmat design_mat_env)
{
    sp_dmat design_mat{std::move(design_mat_env)};
    env_cov_mats_.emplace_back(design_mat * design_mat.t());
    sigma_names_.emplace_back(std::move(name));
}

void GBLUP::set_model()
{
    n_group_effects_ = env_cov_mats_.size();
    n_individuals_ = phenotype_.n_elem;
    n_genetic_effects_ = genetic_cov_mats_.size();
    n_random_effects_ = n_group_effects_ + n_genetic_effects_ + 1;
    sigma_ = arma::zeros(n_group_effects_ + n_genetic_effects_ + 1);
    sigma_names_.emplace_back("e");
}
void GBLUP::reset()
{
    env_cov_mats_.clear();
    genetic_cov_mats_.clear();
    sigma_names_.clear();
}

void GBLUP::add_genetic_effect(std::string name, dmat genetic_covar_mat)
{
    genetic_cov_mats_.emplace_back(std::move(genetic_covar_mat));
    sigma_names_.emplace_back(std::move(name));
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
