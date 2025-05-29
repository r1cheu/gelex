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

void GBLUP::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_mat)
{
    auto n_levels = levels.size();
    fixed_ = FixedEffect{
        std::move(names),
        std::move(levels),
        std::move(design_mat),
        arma::zeros(n_levels)};
}

void GBLUP::add_random_effect(std::string&& name, sp_dmat&& design_mat)
{
    sp_dmat cov = design_mat * design_mat.t();
    random_.add(
        std::move(name),
        effect_type::random,
        std::move(design_mat),
        std::move(cov));
}

void GBLUP::add_genetic_effect(
    std::string&& name,
    sp_dmat&& design_mat,
    const dmat& cov_mat)
{
    if (check_eye(design_mat))
    {
        random_.add(
            std::move(name),
            effect_type::genetic,
            std::move(design_mat),
            cov_mat);
    }
    else
    {
        dmat _cov_mat = design_mat * cov_mat * design_mat.t();
        random_.add(
            std::move(name),
            effect_type::genetic,
            std::move(design_mat),
            std::move(_cov_mat));
    }
}

void GBLUP::add_gxe_effect(
    string&& name,
    sp_dmat&& design_mat_genetic,
    const dmat& genetic_cov_mat,
    const dmat& design_mat)

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

    dmat r_cov = design_mat * design_mat.t();
    g_cov = g_cov % r_cov;
    g_cov = g_cov / arma::trace(g_cov) * static_cast<double>(n_individuals_);

    random_.add(
        std::move(name),
        effect_type::gxe,
        design_mat_genetic,
        std::move(g_cov));
}

void GBLUP::add_residual()
{
    random_.add(
        "e",
        effect_type::residual,
        arma::dmat(),
        arma::speye(n_individuals_, n_individuals_));
}

void GBLUP::clear()
{
    formula_ = "";
    n_individuals_ = 0;
    phenotype_.clear();

    random_.clear();
    fixed_.clear();
}

}  // namespace gelex
