// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_GENOTYPE_GEBV_H_
#define GELEX_BAYES_GENOTYPE_GEBV_H_

#include <Eigen/Core>

#include "gelex/bayes/genotype/projection.h"
#include "gelex/io/csc_reader.h"

namespace gelex
{

// Genomic estimated breeding values of one draw: overwrites `target`
// (n_individuals) with the projection applied to column `draw` of
// `coefficients` (n_markers, n_draws). Requires coefficients.rows() ==
// projection.cols() and target.size() == projection.rows().
auto gebv_draw(
    const bayes::GeneticProjection& projection,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    Eigen::Index draw,
    Eigen::Ref<Eigen::VectorXd> target) -> void;

auto gebv_draw(
    const bayes::GeneticProjection& projection,
    const CscReader::sparse_map_type<double>& coefficients,
    Eigen::Index draw,
    Eigen::Ref<Eigen::VectorXd> target) -> void;

}  // namespace gelex

#endif  // GELEX_BAYES_GENOTYPE_GEBV_H_
