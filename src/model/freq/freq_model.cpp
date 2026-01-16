#include "gelex/model/freq/model.h"
#include "utils/math_utils.h"

namespace gelex
{

FreqModel::FreqModel(DataPipe& data_pipe)
    : phenotype_(std::move(data_pipe).take_phenotype()),
      phenotype_variance_(detail::var(phenotype_)[0]),
      fixed_(std::move(data_pipe).take_fixed_effects())
{
    num_individuals_ = phenotype_.size();

    if (data_pipe.has_additive_grm())
    {
        add_genetic(
            freq::GeneticEffect{
                .name = "Additive",
                .K = std::move(data_pipe).take_additive_grm()});
    }

    if (data_pipe.has_dominance_grm())
    {
        add_genetic(
            freq::GeneticEffect{
                .name = "Dominance",
                .K = std::move(data_pipe).take_dominance_grm()});
    }
}

auto FreqModel::add_random(freq::RandomEffect&& effect) -> void
{
    random_.push_back(std::move(effect));
    effects_view_.emplace_back(&random_.back());
}

auto FreqModel::add_genetic(freq::GeneticEffect&& effect) -> void
{
    genetic_.push_back(std::move(effect));
    effects_view_.emplace_back(&genetic_.back());
}

FreqState::FreqState(const FreqModel& model)
    : phenotype_variance_(model.phenotype_variance()), fixed_(model.fixed())
{
    for (const auto& r : model.random())
    {
        random_.emplace_back(r);
    }
    for (const auto& g : model.genetic())
    {
        genetic_.emplace_back(g);
    }
    init_variance_components(model);
}

auto FreqState::compute_heritability() -> void
{
    double total_genetic_variance = 0.0;
    for (const auto& g : genetic_)
    {
        total_genetic_variance += g.variance;
    }

    double total_random_variance = 0.0;
    for (const auto& r : random_)
    {
        total_random_variance += r.variance;
    }

    phenotype_variance_
        = total_genetic_variance + total_random_variance + residual_.variance;

    if (phenotype_variance_ > 0.0)
    {
        for (auto& g : genetic_)
        {
            g.heritability = g.variance / phenotype_variance_;
        }
    }
}

auto FreqState::init_variance_components(const FreqModel& model) -> void
{
    const double heritability = 0.5;
    const double random_proportion = 0.2;  // exclude genetic random effects
    const double init_residual_variance
        = phenotype_variance_ * (1.0 - heritability);

    double init_genetic_variance = phenotype_variance_ * heritability;
    double init_random_variance = phenotype_variance_ * random_proportion;

    const auto num_genetic = static_cast<double>(model.genetic().size());
    const auto num_random = static_cast<double>(model.random().size());

    init_genetic_variance
        = num_genetic < 1 ? 0.0 : init_genetic_variance / num_genetic;
    init_random_variance
        = num_random < 1 ? 0.0 : init_random_variance / num_random;

    for (auto& g : genetic_)
    {
        g.variance = init_genetic_variance;
    }
    for (auto& r : random_)
    {
        r.variance = init_random_variance;
    }
    residual_.variance = init_residual_variance;
}

}  // namespace gelex
