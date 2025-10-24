#include "gelex/model/bayes/model.h"

#include <optional>
#include <string>

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <Eigen/Core>

#include "../src/data/math_utils.h"
#include "../src/model/bayes/bayes_effects.h"
#include "gelex/data/data_pipe.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/model/effects.h"

namespace gelex
{
using Eigen::Index;

using Eigen::MatrixXd;
using Eigen::VectorXd;

BayesModel::BayesModel(VectorXd&& phenotype, BayesAlphabet alphabet)
    : phenotype_(std::move(phenotype)), alphabet_(alphabet)
{
    num_individuals_ = phenotype_.rows();         // NOLINT
    phenotype_var_ = detail::var(phenotype_)(0);  // NOLINT
}

auto BayesModel::create(DataPipe& data_pipe, BayesAlphabet alphabet)
    -> std::expected<BayesModel, Error>
{
    BayesModel model(std::move(data_pipe).take_phenotype(), alphabet);

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
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
{
    fixed_.emplace(std::move(levels), std::move(design_matrix));
}

void BayesModel::add_random_effect(
    std::vector<std::string>&& levels,
    MatrixXd&& design_matrix)
{
    random_.emplace_back(std::move(levels), std::move(design_matrix));
}

void BayesModel::add_additive(GenotypeMap&& matrix)
{
    additive_.emplace(std::move(matrix));
}

void BayesModel::add_additive(GenotypeMatrix&& matrix)
{
    additive_.emplace(std::move(matrix));
}

void BayesModel::add_dominance(GenotypeMap&& matrix)
{
    dominant_.emplace(std::move(matrix));
}

void BayesModel::add_dominance(GenotypeMatrix&& matrix)
{
    dominant_.emplace(std::move(matrix));
}

BayesState::BayesState(const BayesModel& model)
{
    if (const auto* effect = model.fixed(); effect)
    {
        fixed_.emplace(*effect);
    }

    if (const auto* effect = model.additive(); effect)
    {
        additive_.emplace(*effect, model.is_mixture_model());
    }

    if (const auto* effect = model.dominant(); effect)
    {
        dominant_.emplace(*effect);
    }

    if (const auto& effects = model.random(); !effects.empty())
    {
        random_.reserve(effects.size());
        for (const auto& effect : effects)
        {
            random_.emplace_back(effect);
        }
    }
    residual_.y_adj = model.phenotype().array() - model.phenotype().mean();
    residual_.variance = model.residual().init_variance;
}

void BayesState::compute_heritability()
{
    double sum_var = 0;

    for (const auto& rand : random_)
    {
        sum_var += rand.variance;
    }

    if (additive_)
    {
        sum_var += additive_->variance;
    }
    if (dominant_)
    {
        sum_var += dominant_->variance;
    }

    sum_var += residual_.variance;

    if (additive_)
    {
        additive_->heritability = additive_->variance / sum_var;
    }
    if (dominant_)
    {
        dominant_->heritability = dominant_->variance / sum_var;
    }
}

}  // namespace gelex
