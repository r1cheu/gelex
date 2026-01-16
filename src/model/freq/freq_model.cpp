#include "gelex/model/freq/model.h"

namespace gelex
{

FreqModel::FreqModel(DataPipe& data_pipe)
    : phenotype_(std::move(data_pipe).take_phenotype()),
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

FreqState::FreqState(const FreqModel& model) : fixed_(model.fixed())
{
    for (const auto& r : model.random())
    {
        random_.emplace_back(r);
    }
    for (const auto& g : model.genetic())
    {
        genetic_.emplace_back(g);
    }
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

    phenotypic_variance_
        = total_genetic_variance + total_random_variance + residual_.variance;

    if (phenotypic_variance_ > 0.0)
    {
        for (auto& g : genetic_)
        {
            g.heritability = g.variance / phenotypic_variance_;
        }
    }
}

}  // namespace gelex
