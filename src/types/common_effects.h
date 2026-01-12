#ifndef GELEX_INTERNAL_TYPES_COMMON_EFFECTS_H_
#define GELEX_INTERNAL_TYPES_COMMON_EFFECTS_H_

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct QuantitativeCovariate
{
    std::vector<std::string> names;
    Eigen::MatrixXd X;
};

struct DiscreteCovariate
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    Eigen::MatrixXd X;
};

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

    auto operator[](size_t i) -> CovariateInfoView
    {
        return CovariateInfoView{
            .name = names[i],
            .levels = levels[i] ? std::span<const std::string>(
                                      levels[i]->data(), levels[i]->size())
                                : std::span<const std::string>{},
            .reference_level
            = reference_levels[i] ? *reference_levels[i] : std::string_view{}};
    }

    static auto build_freq(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedEffect;

    static auto build_freq(Eigen::Index n_samples) -> FixedEffect;

    static auto build_bayes(
        std::optional<QuantitativeCovariate> qcovariate,
        std::optional<DiscreteCovariate> dcovariate) -> FixedEffect
    {
        auto fe = build_freq(std::move(qcovariate), std::move(dcovariate));
        fe.cols_norm = fe.X.colwise().squaredNorm();
        return fe;
    }

    static auto build_bayes(Eigen::Index n_samples) -> FixedEffect
    {
        auto fe = build_freq(n_samples);
        fe.cols_norm = fe.X.colwise().squaredNorm();
        return fe;
    }
};

}  // namespace gelex

#endif  // GELEX_INTERNAL_TYPES_COMMON_EFFECTS_H_
