#include "gelex/model/bayes.h"
#include "armadillo"
#include "gelex/model/bayes_policy.h"

#include "gelex/model/effects/bayes_effects.h"

namespace gelex
{
Bayes::Bayes(std::string formula, dvec&& phenotype)
    : formula_(std::move(formula)), phenotype_(std::move(phenotype))
{
    n_individuals_ = phenotype_.n_elem;      // NOLINT
    mu_.value = arma::mean(phenotype_);      // NOLINT
    phenotype_var_ = arma::var(phenotype_);  // NOLINT
    residual_.value = 0.01;
}

void Bayes::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    dmat&& design_mat)
{
    fixed_ = bayes::FixedEffect(
        std::move(names), std::move(levels), std::move(design_mat));
}

void Bayes::add_random_effect(std::string&& name, dmat&& design_mat)
{
    random_.add(
        bayes::RandomEffect(std::move(name), dvec{0.0}, std::move(design_mat)));
}

void Bayes::add_genetic_effect(
    std::string&& name,
    dmat&& genotype,
    BayesAlphabet type)
{
    auto sigma = bayes_trait.at(to_index(type))(genotype);
    auto [pi, num] = bayes_trait_pi.at(to_index(type))();

    genetic_.add(
        bayes::GeneticEffect(
            std::move(name),
            std::move(sigma),
            std::move(genotype),
            bayes::BayesParam(type, std::move(pi), std::move(num))));
}

}  // namespace gelex
