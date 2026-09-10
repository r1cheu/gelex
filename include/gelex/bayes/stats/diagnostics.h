// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_STATS_DIAGNOSTICS_H_
#define GELEX_BAYES_STATS_DIAGNOSTICS_H_

#include <Eigen/Core>

namespace gelex
{

// Read-only (n_draws, n_chains) view with arbitrary strides: each column is
// one chain. Binds a VectorXd, mat.col(i), mat.row(i).transpose(), a block or
// a transposed matrix without copying.
using ChainsView = Eigen::Ref<
    const Eigen::MatrixXd,
    0,
    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>;

struct ChainDiagnostics
{
    double mean;
    double sd;
    double median;
    double hpdi_lower;
    double hpdi_upper;
    double ess;
    double mcse;
    double split_rhat;
};

/**
 * @brief Computes every diagnostic of one parameter from its chains, each
 * stored as a column: pooled sample mean and standard deviation, median,
 * HPDI bounds (narrowest interval with probability mass `prob`), effective
 * sample size (Geyer initial monotone sequence), Monte Carlo standard error,
 * and split R-hat over the halves of every chain. Requires at least 4 draws
 * per chain; a chain that never moves gets NaN for ess, mcse and split_rhat.
 *
 * @param draws (n_draws, n_chains); left untouched
 * @param prob probability mass of the HPDI, in (0, 1]
 */
auto diagnose_chain(const ChainsView& draws, double prob = 0.95)
    -> ChainDiagnostics;

template <typename Derived>
auto diagnose_chain(const Eigen::DenseBase<Derived>& draws, double prob = 0.95)
    -> ChainDiagnostics
{
    if constexpr (
        Derived::IsVectorAtCompileTime && Derived::RowsAtCompileTime == 1)
    {
        return diagnose_chain(ChainsView{draws.derived().transpose()}, prob);
    }
    else
    {
        return diagnose_chain(ChainsView{draws.derived()}, prob);
    }
}

}  // namespace gelex

#endif  // GELEX_BAYES_STATS_DIAGNOSTICS_H_
