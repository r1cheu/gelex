#pragma once
#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace gelex::freq
{

struct QCovariateEffect
{
    std::vector<std::string> names;
    Eigen::MatrixXd X;
};

struct DCovariateEffect
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    Eigen::MatrixXd X;
};

struct CovariateInfoView
{
    std::string_view name;
    std::span<const std::string> levels;
    std::string_view reference_level;
};

struct FixedEffect
{
    std::vector<std::string> names;
    std::vector<std::optional<std::vector<std::string>>> levels;
    std::vector<std::optional<std::string>> reference_levels;
    Eigen::MatrixXd X;

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

    static FixedEffect build(
        QCovariateEffect&& qcovariate,
        DCovariateEffect&& dcovariate);
};

struct RandomEffect
{
    std::string name;
    std::vector<std::string> levels;
    Eigen::SparseMatrix<double> Z;  // incidence matrix
    Eigen::SparseMatrix<double> K;  // covariance matrix
};

struct GeneticEffect
{
    std::string name;
    Eigen::MatrixXd K;  // GRM matrix
};

struct FixedState
{
    explicit FixedState(const FixedEffect& effect);
    Eigen::VectorXd coeff;
    Eigen::VectorXd se;
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect);
    std::string name;
    Eigen::VectorXd blup;
    double variance;
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect);
    std::string name;
    Eigen::VectorXd ebv;
    double variance;
};

struct ResidualState
{
    double variance{};
    double heritability{};
};

}  // namespace gelex::freq
