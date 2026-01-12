#include "../src/types/fixed_effects.h"

#include <iterator>
#include <optional>

#include "gelex/exception.h"

namespace gelex
{

auto FixedEffect::build(
    std::optional<QuantitativeCovariate> qcovariate,
    std::optional<DiscreteCovariate> dcovariate) -> FixedEffect
{
    FixedEffect fe;

    Eigen::Index n_samples = 0;
    if (qcovariate)
    {
        n_samples = qcovariate->X.rows();
    }
    else if (dcovariate)
    {
        n_samples = dcovariate->X.rows();
    }
    else
    {
        throw gelex::InvalidInputException(
            "At least one covariate must be provided");
    }

    const auto qcov_cols = qcovariate ? qcovariate->X.cols() : 0;
    const auto dcov_cols = dcovariate ? dcovariate->X.cols() : 0;
    const auto n_cols = 1 + qcov_cols + dcov_cols;

    fe.names.reserve(n_cols);
    fe.names.emplace_back("Intercept");

    auto move_insert = [](auto& container, auto begin, auto end)
    {
        container.insert(
            container.end(),
            std::make_move_iterator(begin),
            std::make_move_iterator(end));
    };

    if (qcovariate)
    {
        move_insert(
            fe.names, qcovariate->names.begin(), qcovariate->names.end());
    }
    if (dcovariate)
    {
        move_insert(
            fe.names, dcovariate->names.begin(), dcovariate->names.end());
    }

    fe.levels.reserve(n_cols);
    fe.levels.emplace_back(std::nullopt);
    if (qcovariate)
    {
        fe.levels.insert(
            fe.levels.end(), qcovariate->names.size(), std::nullopt);
    }
    if (dcovariate)
    {
        move_insert(
            fe.levels, dcovariate->levels.begin(), dcovariate->levels.end());
    }

    fe.reference_levels.reserve(n_cols);
    fe.reference_levels.emplace_back(std::nullopt);
    if (qcovariate)
    {
        fe.reference_levels.insert(
            fe.reference_levels.end(), qcovariate->names.size(), std::nullopt);
    }
    if (dcovariate)
    {
        move_insert(
            fe.reference_levels,
            dcovariate->reference_levels.begin(),
            dcovariate->reference_levels.end());
    }

    fe.X = Eigen::MatrixXd::Zero(n_samples, n_cols);
    fe.X.col(0).setOnes();

    Eigen::Index col_offset = 1;
    if (qcovariate)
    {
        fe.X.middleCols(col_offset, qcov_cols) = qcovariate->X;
        col_offset += qcov_cols;
    }
    if (dcovariate)
    {
        fe.X.middleCols(col_offset, dcov_cols) = dcovariate->X;
    }

    return fe;
}

auto FixedEffect::build(Eigen::Index n_samples) -> FixedEffect
{
    FixedEffect fe;

    const Eigen::Index n_cols = 1;

    fe.names.reserve(n_cols);
    fe.names.emplace_back("Intercept");

    fe.levels.reserve(n_cols);
    fe.levels.emplace_back(std::nullopt);

    fe.reference_levels.reserve(n_cols);
    fe.reference_levels.emplace_back(std::nullopt);

    fe.X = Eigen::MatrixXd::Zero(n_samples, n_cols);
    fe.X.col(0).setOnes();

    return fe;
}

}  // namespace gelex
