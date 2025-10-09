#include "gelex/model/bayes/model.h"

#include <memory>
#include <random>
#include <stdexcept>
#include <string>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/math_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/utils/formatter.h"
#include "model/bayes/distribution.h"
#include "model/bayes/traits/base_trait.h"

namespace gelex
{
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

BayesModel::BayesModel(VectorXd&& phenotype, BayesAlphabet type)
    : phenotype_(std::move(phenotype))
{
    n_individuals_ = phenotype_.rows();           // NOLINT
    phenotype_var_ = detail::var(phenotype_)(0);  // NOLINT
    model_trait_ = create_genetic_trait(type);
    set_sigma_prior_manual("e", 4, 0, phenotype_var_ * 0.5);  // NOLINT
}

auto BayesModel::create_from_datapipe(DataPipe& data_pipe, BayesAlphabet type)
    -> std::expected<BayesModel, Error>
{
    // Create model with phenotype from DataPipe
    BayesModel model(std::move(data_pipe).take_phenotype(), type);

    // Add fixed effects if available
    if (data_pipe.has_fixed_effects())
    {
        auto fixed_effect_names = data_pipe.fixed_effect_names();
        model.add_fixed_effect(
            std::move(fixed_effect_names),
            std::move(data_pipe).take_fixed_effects());
    }

    return model;
}

void BayesModel::add_fixed_effect(
    std::vector<std::string>&& names,
    MatrixXd&& design_matrix)
{
    fixed_ = std::make_unique<bayes::FixedEffect>(
        std::move(names), std::move(design_matrix));
}

void BayesModel::add_random_effect(std::string&& name, MatrixXd&& design_matrix)
{
    random_.add(std::move(name), std::move(design_matrix));
    for (auto& effect : random_.effects())
    {
        set_sigma_prior(effect.name, 4, 0.5);
    }
    set_sigma_prior_manual("e", 4, 0, phenotype_var_ * 0.5);
}

void BayesModel::add_additive_effect(GenotypeMap&& matrix)

{
    Index n_snp = matrix.cols();
    additive_ = std::make_unique<bayes::AdditiveEffect>(
        "additive",
        std::move(matrix),
        model_trait_->default_marker_variance(n_snp),
        model_trait_->default_pi());
    set_sigma_prior(additive_->name, 4.0, 0.5);
}

void BayesModel::add_dominance_effect(GenotypeMap&& matrix)
{
    Index n_snp = matrix.cols();
    dominant_ = std::make_unique<bayes::DominantEffect>(
        "dominant", std::move(matrix), 0, 2);
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
        random_effect->effect_variance.setConstant(init_sigma);
        return;
    }
    if (additive_ && additive_->name == name)
    {
        additive_->prior = {nu, s2};
        additive_->marker_variance.setConstant(init_sigma);
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
        init_sigma = var / additive_->design_matrix.variance().sum()
                     / (1 - additive_->pi(0));
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
        sum_var += rand.effect_variance(0);
    }

    if (additive)
    {
        sum_var += additive->effect_variance;
    }
    if (dominant)
    {
        sum_var += dominant->effect_variance;
    }
    sum_var += residual.value;

    if (additive)
    {
        additive->heritability = additive->effect_variance / sum_var;
    }
    if (dominant)
    {
        dominant->heritability = dominant->effect_variance / sum_var;
    }
}

}  // namespace gelex
