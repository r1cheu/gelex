#ifndef GELEX_MODEL_FREQ_FREQMODEL_H_
#define GELEX_MODEL_FREQ_FREQMODEL_H_

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "../src/types/fixed_effects.h"
#include "../src/types/freq_effect.h"
#include "gelex/data/data_pipe.h"

namespace gelex
{

class FreqModel
{
   public:
    explicit FreqModel(DataPipe& data_pipe);
    auto fixed() const -> const FixedEffect& { return fixed_; }
    auto fixed() -> FixedEffect& { return fixed_; }

    auto random() const -> const std::vector<freq::RandomEffect>&
    {
        return random_;
    }
    auto random() -> std::vector<freq::RandomEffect>& { return random_; }

    auto genetic() const -> const std::vector<freq::GeneticEffect>&
    {
        return genetic_;
    }
    auto genetic() -> std::vector<freq::GeneticEffect>& { return genetic_; }

    // return random and genetic effects as a single vector of pointers, for
    // convenience
    auto effects_view() const -> const
        std::vector<std::variant<freq::RandomEffect*, freq::GeneticEffect*>>&
    {
        return effects_view_;
    }

    auto phenotype() const -> const Eigen::VectorXd& { return phenotype_; }
    auto phenotype_variance() const -> double { return phenotype_variance_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

   private:
    auto add_random(freq::RandomEffect&& effect) -> void;
    auto add_genetic(freq::GeneticEffect&& effect) -> void;

    Eigen::Index num_individuals_{};

    Eigen::VectorXd phenotype_;
    double phenotype_variance_;

    FixedEffect fixed_;
    std::vector<freq::RandomEffect> random_;
    std::vector<freq::GeneticEffect> genetic_;

    std::vector<std::variant<freq::RandomEffect*, freq::GeneticEffect*>>
        effects_view_;
};

class FreqState
{
   public:
    explicit FreqState(const FreqModel& model);

    freq::FixedState& fixed() { return fixed_; }
    const freq::FixedState& fixed() const { return fixed_; }

    std::vector<freq::RandomState>& random() { return random_; }
    const std::vector<freq::RandomState>& random() const { return random_; }

    std::vector<freq::GeneticState>& genetic() { return genetic_; }
    const std::vector<freq::GeneticState>& genetic() const { return genetic_; }

    freq::ResidualState& residual() { return residual_; }
    const freq::ResidualState& residual() const { return residual_; }

    void compute_heritability();

   private:
    double phenotype_variance_;
    freq::FixedState fixed_;
    std::vector<freq::RandomState> random_;
    std::vector<freq::GeneticState> genetic_;
    freq::ResidualState residual_;

    void init_variance_components(const FreqModel& model);
};

}  // namespace gelex

#endif  // GELEX_MODEL_FREQ_FREQMODEL_H_
