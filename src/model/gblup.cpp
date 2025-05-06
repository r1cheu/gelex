#include "gelex/model/gblup.h"

#include <fmt/format.h>
#include <armadillo>
#include <utility>

#include "gelex/utils.h"

namespace gelex
{

GBLUP::GBLUP(std::string formula, dvec phenotype)
    : formula_{std::move(formula)}, phenotype_{std::move(phenotype)}
{
    n_individuals_ = phenotype_.n_elem;
}

void GBLUP::add_group_effect(std::string name, sp_dmat design_mat_group)
{
    sp_dmat cov = design_mat_group * design_mat_group.t();
    effects_.add_effect(
        std::move(name),
        effect_type::group,
        std::move(design_mat_group),
        std::move(cov));
}
void GBLUP::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_mat_beta)
{
    design_mat_beta_ = std::move(design_mat_beta);
    fixed_eff_names_ = std::move(names);
    fixed_eff_levels_ = std::move(levels);
    n_common_effects_ = design_mat_beta_.n_cols;
    beta_ = arma::zeros(n_common_effects_);
}

void GBLUP::clear()
{
    effects_.clear();
}

void GBLUP::add_genetic_effect(
    std::string name,
    sp_dmat design_mat_genetic,
    const dmat& genetic_covar_mat)
{
    if (check_eye(design_mat_genetic))
    {
        effects_.add_effect(
            std::move(name),
            effect_type::genetic,
            std::move(design_mat_genetic),
            genetic_covar_mat);
    }
    else
    {
        dmat cov_mat
            = design_mat_genetic * genetic_covar_mat * design_mat_genetic.t();
        effects_.add_effect(
            std::move(name),
            effect_type::genetic,
            std::move(design_mat_genetic),
            std::move(cov_mat));
    }
}

void GBLUP::add_gxe_effect(
    string name,
    sp_dmat design_mat_genetic,
    const dmat& genetic_cov_mat,
    const dmat& design_mat_group)

{
    sp_dmat g_design = std::move(design_mat_genetic);
    dmat g_cov;
    if (check_eye(g_design))
    {
        g_cov = genetic_cov_mat;
    }
    else
    {
        g_cov = g_design * genetic_cov_mat * g_design.t();
    }

    dmat r_cov = design_mat_group * design_mat_group.t();
    g_cov = g_cov % r_cov;
    g_cov = g_cov / arma::trace(g_cov) * static_cast<double>(n_individuals_);

    effects_.add_effect(
        std::move(name),
        effect_type::gxe,
        design_mat_genetic,
        std::move(g_cov));
}

}  // namespace gelex
