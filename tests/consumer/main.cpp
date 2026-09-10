// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <Eigen/Core>
#include <gelex/bayes/stats/diagnostics.h>
#include <gelex/data/encode/matrix.h>
#include <gelex/data/genotype_method.h>
#include <gelex/genetic_mode.h>

auto main() -> int
{
    const auto modes = gelex::GeneticMode::A | gelex::GeneticMode::D;
    Eigen::MatrixXd genotypes{{0.0}, {1.0}, {2.0}};
    gelex::encode_inplace(
        genotypes, gelex::GeneticMode::A, gelex::GenotypeMethod::Center);
    const Eigen::MatrixXd expected{{-1.0}, {0.0}, {1.0}};

    const Eigen::VectorXd draws = Eigen::VectorXd::LinSpaced(8, 0.0, 7.0);
    const auto diagnostics = gelex::diagnose_chain(draws);

    return modes.size() == 2 && diagnostics.mean == 3.5
                   && genotypes.isApprox(expected)
               ? 0
               : 1;
}
