#include "gelex/model/freq/model.h"

#include <armadillo>

namespace gelex
{

using arma::dmat;
using arma::dvec;
using arma::sp_dmat;

GBLUP::GBLUP(std::string formula, dvec&& phenotype)
    : formula_{std::move(formula)}, phenotype_{std::move(phenotype)}
{
    n_individuals_ = phenotype_.n_elem;
    effects_.add("e", effect_type::residual, 0, 0);
}

void GBLUP::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_matrix)
{
    fixed_ = std::make_unique<FixedEffect>(
        std::move(names), std::move(levels), std::move(design_matrix));
}

void GBLUP::add_random_effect(std::string&& name, sp_dmat&& design_matrix)
{
    random_.add(std::move(name), std::move(design_matrix));
    effects_.add(random_.back().name, effect_type::random, 0, 0);
}

void GBLUP::add_genetic_effect(
    std::string&& name,
    sp_dmat&& design_matrix,
    const dmat& genetic_relationship_matrix)
{
    genetic_.add(
        std::move(name), std::move(design_matrix), genetic_relationship_matrix);
    effects_.add(genetic_.back().name, effect_type::genetic, 0, 0);
}

void GBLUP::add_gxe_effect(
    std::string&& name,
    sp_dmat&& genetic_design_matrix,
    const dmat& genetic_relationship_matrix,
    arma::sp_dmat&& env_design_matrix)
{
    gxe_.add(
        std::move(name),
        std::move(genetic_design_matrix),
        genetic_relationship_matrix,
        std::move(env_design_matrix));
    effects_.add(gxe_.back().name, effect_type::gxe, 0, 0);
}

void GBLUP::clear()
{
    formula_ = "";
    n_individuals_ = 0;
    phenotype_.clear();

    fixed_.reset();
    random_.clear();
    genetic_.clear();
    gxe_.clear();
    effects_.clear();
    residual_.sigma = 0;
}

dvec GBLUP::var_comp() const
{
    std::vector<double> var_comp_values;
    var_comp_values.emplace_back(residual_.sigma);
    auto add_var_comp = [&](const auto& effects)
    {
        for (const auto& eff : effects)
        {
            var_comp_values.emplace_back(eff.sigma);
        }
    };
    add_var_comp(random_);
    add_var_comp(genetic_);
    add_var_comp(gxe_);
    return dvec(var_comp_values);
}

size_t GBLUP::n_var_comp() const
{
    return random_.size() + genetic_.size() + gxe_.size() + 1;
}

void GBLUP::init_var_comp(double var)
{
    var /= static_cast<double>(n_var_comp());
    arma::dvec var_comp(n_var_comp(), arma::fill::value(var));
    set_var_comp(var_comp);
}

void GBLUP::set_var_comp(const dvec& var_comp)
{
    size_t idx{0};
    residual_.sigma = var_comp.at(idx++);
    effects_.get(residual_.name)->sigma = residual_.sigma;

    auto set_sigma = [&](auto& effects)
    {
        for (auto& eff : effects)
        {
            double sigma = var_comp.at(idx++);
            eff.sigma = sigma;
            effects_.get(eff.name)->sigma = sigma;
        }
    };
    set_sigma(random_);
    set_sigma(genetic_);
    set_sigma(gxe_);
}
}  // namespace gelex
