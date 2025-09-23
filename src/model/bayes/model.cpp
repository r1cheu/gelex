#include "gelex/model/bayes/model.h"

#include <memory>
#include <stdexcept>
#include <string>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/genotype_mmap.h"
#include "../src/data/math_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "gelex/utils/formatter.h"

namespace gelex
{
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

BayesModel::BayesModel(VectorXd&& phenotype) : phenotype_(std::move(phenotype))
{
    n_individuals_ = phenotype_.rows();                        // NOLINT
    phenotype_var_ = detail::var(phenotype_)(0);               // NOLINT
    set_sigma_prior_manual("e", -2, 0, phenotype_var_ * 0.5);  // NOLINT
}

void BayesModel::add_fixed_effect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
{
    fixed_ = std::make_unique<bayes::FixedEffect>(
        std::move(names), std::move(levels), std::move(design_matrix));
}

void BayesModel::add_random_effect(std::string&& name, MatrixXd&& design_matrix)
{
    random_.add(std::move(name), std::move(design_matrix));
    for (auto& effect : random_.effects())
    {
        set_sigma_prior(effect.name, 4, 0.5);
    }
    set_sigma_prior_manual("e", -2, 0, phenotype_var_ * 0.5);
}

void BayesModel::add_additive_effect(
    detail::GenotypeMap&& matrix,
    BayesAlphabet type)

{
    Index n_snp = matrix.cols();
    additive_ = std::make_unique<bayes::AdditiveEffect>(
        "add",
        std::move(matrix),
        type,
        model_trait_->default_sigma(n_snp),
        model_trait_->default_pi());
    set_sigma_prior(additive_->name, 4.0, 0.5);
}

void BayesModel::add_dominance_effect(
    detail::GenotypeMap&& matrix,
    BayesAlphabet type)
{
    Index n_snp = matrix.cols();
    dominant_ = std::make_unique<bayes::DominantEffect>(
        "dom", std::move(matrix), 0, 2);
    set_sigma_prior(dominant_->name, 4.0, 0.5);
}

void BayesModel::set_sigma_prior_manual(
    const std::string& name,
    double nu,
    double s2,
    double init_sigma)
{
    if (auto* random_effect = random_.get(name))
    {
        random_effect->prior = {nu, s2};
        random_effect->sigma.setConstant(init_sigma);
        return;
    }
    if (additive_ && additive_->name == name)
    {
        additive_->prior = {nu, s2};
        additive_->sigma.setConstant(init_sigma);
        return;
    }
    if (name == "e")
    {
        residual_.prior = {nu, s2};
        residual_.value = init_sigma;
        return;
    }

    // Build available effects list safely
    std::vector<std::string> available_effects;
    for (const auto& effect_name : random_.keys())
    {
        available_effects.push_back(effect_name);
    }
    if (additive_)
    {
        available_effects.push_back(additive_->name);
    }
    available_effects.emplace_back("e");

    throw std::runtime_error(
        fmt::format(
            "Effect not found: '{}'. Available effects: {}",
            name,
            fmt::join(available_effects, ", ")));
}

void BayesModel::set_sigma_prior(
    const std::string& name,
    double nu,
    double prop)
{
    double var = prop * phenotype_var_;
    size_t nr = 0;
    double init_sigma = 0.0;
    double s2 = 0.0;

    for (const auto& effect : random_.effects())
    {
        nr += effect.design_matrix.cols();
    }

    if (additive_ && additive_->name == name)
    {
        auto n_snp = static_cast<double>(additive_->design_matrix.cols());
        auto p1 = additive_->pi(1);
        init_sigma = var / n_snp / p1;
        s2 = (nu - 2) / nu * init_sigma;
    }
    else
    {
        init_sigma = var / static_cast<double>(nr + 1);
        s2 = (nu - 2) / nu * init_sigma;
    }
    set_sigma_prior_manual(name, nu, s2, init_sigma);
}

// New implementation with default calculation
void BayesModel::set_pi_prior(const std::string& name, const VectorXd& pi)
{
    if (additive_->name == name)
    {
        additive_->pi = pi;
        return;
    }
    throw std::runtime_error(
        fmt::format(
            "Additive effect not found: {}, `{}` are available.",
            name,
            additive_->name));
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

    if (additive_)
    {
        auto prior_str = model_trait_->prior_info(
            additive_->prior.nu, additive_->prior.s2, additive_->pi);
        for (size_t i{}; i < prior_str.size(); ++i)
        {
            if (i == 0)
            {
                summary += fmt::format("{}: {}", additive_->name, prior_str[i]);
            }
            else
            {
                summary += prior_str[i];
            }
            summary += '\n';
        }
    }

    summary += fmt::format(
        "e: {}", sigma_prior("", residual_.prior.nu, residual_.prior.s2));
    summary += '\n';

    return summary;
}

void BayesStatus::compute_heritability()
{
    double sum_var = 0;

    for (const auto& rand : random)
    {
        sum_var += rand.sigma(0);
    }

    if (additive)
    {
        sum_var += additive->variance;
    }
    if (dominant)
    {
        sum_var += dominant->variance;
    }
    sum_var += residual.value;

    if (additive)
    {
        additive->heritability = additive->variance / sum_var;
    }
    if (dominant)
    {
        dominant->heritability = dominant->variance / sum_var;
    }
}

}  // namespace gelex
