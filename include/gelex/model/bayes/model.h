#ifndef GELEX_MODEL_BAYES_MODEL_H_
#define GELEX_MODEL_BAYES_MODEL_H_

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/types/bayes_effects.h"
#include "gelex/data/data_pipe.h"

namespace gelex
{

class GenotypeMap;
class GenotypeMatrix;

class BayesModel
{
   public:
    explicit BayesModel(DataPipe& data_pipe);

    const bayes::FixedEffect* fixed() const
    {
        return fixed_.has_value() ? &fixed_.value() : nullptr;
    }

    bayes::FixedEffect* fixed()
    {
        return fixed_.has_value() ? &fixed_.value() : nullptr;
    }

    const std::vector<bayes::RandomEffect>& random() const { return random_; }
    std::vector<bayes::RandomEffect>& random() { return random_; }

    const bayes::AdditiveEffect* additive() const
    {
        return additive_.has_value() ? &additive_.value() : nullptr;
    }
    bayes::AdditiveEffect* additive()
    {
        return additive_.has_value() ? &additive_.value() : nullptr;
    }

    const bayes::DominantEffect* dominant() const
    {
        return dominant_.has_value() ? &dominant_.value() : nullptr;
    }
    bayes::DominantEffect* dominant()
    {
        return dominant_.has_value() ? &dominant_.value() : nullptr;
    }

    const bayes::Residual& residual() const { return residual_; }
    bayes::Residual& residual() { return residual_; }

    const Eigen::VectorXd& phenotype() const { return phenotype_; }

    double phenotype_variance() const { return phenotype_var_; }
    Eigen::Index num_individuals() const { return num_individuals_; }

   private:
    void add_additive(GenotypeMap&& matrix);
    void add_additive(GenotypeMatrix&& matrix);
    void add_dominance(GenotypeMap&& matrix);
    void add_dominance(GenotypeMatrix&& matrix);

    void add_fixed_effect(
        std::vector<std::string>&& levels,
        Eigen::MatrixXd&& design_matrix);
    void add_random_effect(
        std::vector<std::string>&& levels,
        Eigen::MatrixXd&& design_matrix);

    Eigen::Index num_individuals_{};
    double phenotype_var_{};

    Eigen::VectorXd phenotype_;

    std::optional<bayes::FixedEffect> fixed_;
    std::vector<bayes::RandomEffect> random_;
    std::optional<bayes::AdditiveEffect> additive_;
    std::optional<bayes::DominantEffect> dominant_;
    bayes::Residual residual_;
};

class BayesState
{
   public:
    explicit BayesState(const BayesModel& model);

    bayes::FixedState* fixed() { return fixed_ ? &fixed_.value() : nullptr; }
    const bayes::FixedState* fixed() const
    {
        return fixed_ ? &fixed_.value() : nullptr;
    }
    std::vector<bayes::RandomState>& random() { return random_; }
    const std::vector<bayes::RandomState>& random() const { return random_; }

    bayes::AdditiveState* additive()
    {
        return additive_ ? &additive_.value() : nullptr;
    }
    const bayes::AdditiveState* additive() const
    {
        return additive_ ? &additive_.value() : nullptr;
    }

    bayes::DominantState* dominant()
    {
        return dominant_ ? &dominant_.value() : nullptr;
    }
    const bayes::DominantState* dominant() const
    {
        return dominant_ ? &dominant_.value() : nullptr;
    }

    bayes::ResidualState& residual() { return residual_; }
    const bayes::ResidualState& residual() const { return residual_; }

    void compute_heritability();

   private:
    std::optional<bayes::FixedState> fixed_;
    std::vector<bayes::RandomState> random_;
    std::optional<bayes::AdditiveState> additive_;
    std::optional<bayes::DominantState> dominant_;
    bayes::ResidualState residual_;
};

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_MODEL_H_
