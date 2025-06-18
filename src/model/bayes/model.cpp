#include "gelex/model/bayes/model.h"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <armadillo>

#include "gelex/model/bayes/effects.h"
#include "gelex/model/bayes/policy.h"
#include "gelex/utils/formatter.h"

namespace gelex
{
using arma::dmat;
using arma::dvec;

BayesModel::BayesModel(std::string formula, dvec&& phenotype)
    : formula_(std::move(formula)), phenotype_(std::move(phenotype))
{
    n_individuals_ = phenotype_.n_elem;      // NOLINT
    mu_.value = arma::mean(phenotype_);      // NOLINT
    phenotype_var_ = arma::var(phenotype_);  // NOLINT
    residual_.value = 0.01;
}

void BayesModel::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_mat)
{
    fixed_ = std::make_unique<FixedEffectDesign>(
        std::move(names), std::move(levels), std::move(design_mat));
}

void BayesModel::add_random_effect(std::string&& name, dmat&& design_mat)
{
    random_.add(RandomEffectDesign(std::move(name), std::move(design_mat)));
}

void BayesModel::add_genetic_effect(
    std::string&& name,
    dmat&& genotype,
    BayesAlphabet type)
{
    auto pi = bayes_trait_pi.at(to_index(type))();
    auto sigma = bayes_trait_sigma.at(to_index(type))(genotype);

    genetic_.add(GeneticEffectDesign(
        std::move(name),
        std::move(genotype),
        type,
        std::move(sigma),
        std::move(pi)));
}

void BayesModel::set_sigma_prior(const std::string& name, double nu, double s2)
{
    if (auto* random_effect = random_.get(name))
    {
        random_effect->prior = {nu, s2};
        return;
    }
    if (auto* genetic_effect = genetic_.get(name))
    {
        genetic_effect->prior = {nu, s2};
        return;
    }

    if (name == "e")
    {
        residual_.prior = {nu, s2};
        return;
    }

    throw std::runtime_error(
        fmt::format(
            "Effect not found: {}, `{} {} e` are available.",
            name,
            fmt::join(random_.names(), ", "),
            fmt::join(genetic_.names(), ", ")));
}

void BayesModel::set_pi_prior(const std::string& name, const arma::dvec& pi)
{
    if (auto* effect = genetic_.get(name))
    {
        effect->pi = pi;
        return;
    }
    throw std::runtime_error(
        fmt::format(
            "Genetic effect not found: {}, `{}` are available.",
            name,
            fmt::join(genetic_.names(), ", ")));
}

std::string BayesModel::prior_summary() const
{
    std::string summary{"prior summary:\n"};

    if (random_)
    {
        for (const auto& effect : random_.effects())
        {
            summary += fmt::format(
                "{}: {}",
                effect.name,
                sigma_prior(effect.name, effect.prior.nu, effect.prior.s2));
            summary += '\n';
        }
    }

    if (genetic_)
    {
        for (const auto& effect : genetic_.effects())
        {
            auto prior_str = bayes_trait_prior_str[to_index(effect.type)](
                effect.prior.nu, effect.prior.s2, effect.pi);
            for (size_t i{}; i < prior_str.size(); ++i)
            {
                if (i == 0)
                {
                    summary += fmt::format("{}: {}", effect.name, prior_str[i]);
                }
                else
                {
                    summary += prior_str[i];
                }
                summary += '\n';
            }
        }
    }

    summary += fmt::format(
        "e: {}", sigma_prior("", residual_.prior.nu, residual_.prior.s2));
    summary += '\n';

    return summary;
}

}  // namespace gelex
