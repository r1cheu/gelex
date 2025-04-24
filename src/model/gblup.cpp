#include "gelex/model/gblup.h"

#include <fmt/format.h>
#include <armadillo>
#include <utility>

#include "gelex/utils.h"

namespace gelex
{

GBLUP::GBLUP(std::string formula, dvec phenotype, dmat design_mat_beta)
    : formula_{std::move(formula)},
      phenotype_{std::move(phenotype)},

      design_mat_beta_{std::move(design_mat_beta)}
{
    n_common_effects_ = design_mat_beta_.n_cols;
    beta_ = arma::zeros(n_common_effects_);
    n_individuals_ = phenotype_.n_elem;
}

void GBLUP::add_group_effect(std::string name, dmat design_mat_group)
{
    sp_dmat cov = sp_dmat{design_mat_group * design_mat_group.t()};
    effects_.add_effect(
        std::move(name),
        effect_type::group,
        std::move(design_mat_group),
        std::move(cov));
}

void GBLUP::clear()
{
    effects_.clear();
}

void GBLUP::add_genetic_effect(
    std::string name,
    sp_dmat design_mat_genetic,
    dmat genetic_covar_mat)
{
    if (check_eye(design_mat_genetic))
    {
        effects_.add_effect(
            std::move(name),
            effect_type::genetic,
            std::move(design_mat_genetic),
            std::move(genetic_covar_mat));
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
    const string& genetic_name,
    const dvec& design_mat_group)
{
    auto* genetic_effect = effects_.get(genetic_name);
    if (genetic_effect == nullptr)
    {
        throw std::runtime_error(
            fmt::format(
                "Genetic effect {} not found in GBLUP model", genetic_name));
    }

    const dmat& genetic_cov = std::get<dmat>(genetic_effect->cov_mat);
    const sp_dmat& design_mat = std::get<sp_dmat>(genetic_effect->design_mat);

    dmat cov = diagmat(design_mat_group) * genetic_cov
               * diagmat(design_mat_group).t();
    if (check_eye(design_mat))
    {
        effects_.add_effect(
            std::move(name),
            effect_type::gxe,
            sp_dmat(arma::diagmat(design_mat_group)),
            std::move(cov));
    }
    else
    {
        sp_dmat design_mat_gxe
            = sp_dmat(diagmat(design_mat_group) * design_mat);
        effects_.add_effect(
            std::move(name),
            effect_type::gxe,
            std::move(design_mat_gxe),
            std::move(cov));
    }
}

}  // namespace gelex
