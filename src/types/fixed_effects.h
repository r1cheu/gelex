#ifndef GELEX_TYPES_FIXED_EFFECTS_H_
#define GELEX_TYPES_FIXED_EFFECTS_H_

#include <optional>
#include <span>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "covariates.h"

namespace gelex
{

struct FixedEffect
{
    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    std::vector<std::optional<std::string>> reference_levels;
    Eigen::MatrixXd X;
    std::optional<Eigen::VectorXd> cols_norm;

    struct CovariateInfoView
    {
        std::string_view name;
        std::span<const std::string> levels;
        std::string_view reference_level;
    };

    CovariateInfoView operator[](size_t i)
    {
        return CovariateInfoView{
            .name = names[i],
            .levels = levels[i] ? std::span<const std::string>(
                                      levels[i]->data(), levels[i]->size())
                                : std::span<const std::string>{},
            .reference_level
            = reference_levels[i] ? *reference_levels[i] : std::string_view{}};
    }

    static auto build(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedEffect;

    static auto build(Eigen::Index n_samples) -> FixedEffect;

    static auto build_bayes(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedEffect
    {
        auto fe = build(std::move(qcovariate), std::move(dcovariate));
        fe.cols_norm = fe.X.colwise().squaredNorm();
        return fe;
    }

    static auto build_bayes(Eigen::Index n_samples) -> FixedEffect
    {
        auto fe = build(n_samples);
        fe.cols_norm = fe.X.colwise().squaredNorm();
        return fe;
    }
};

}  // namespace gelex

#endif  // GELEX_TYPES_FIXED_EFFECTS_H_
