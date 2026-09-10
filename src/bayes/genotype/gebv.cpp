// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/bayes/genotype/gebv.h"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cassert>

#include "gelex/bayes/genotype/projection.h"
#include "gelex/io/csc_reader.h"

namespace gelex
{

auto gebv_draw(
    const bayes::GeneticProjection& projection,
    const Eigen::Ref<const Eigen::MatrixXd>& coefficients,
    Eigen::Index draw,
    Eigen::Ref<Eigen::VectorXd> target) -> void
{
    assert(coefficients.rows() == projection.cols());
    assert(target.size() == projection.rows());
    target.setZero();
    for (Eigen::Index marker = 0; marker < coefficients.rows(); ++marker)
    {
        const double coefficient = coefficients(marker, draw);
        if (coefficient != 0.0)
        {
            projection.axpy(marker, coefficient, target);
        }
    }
}

auto gebv_draw(
    const bayes::GeneticProjection& projection,
    const CscReader::sparse_map_type<double>& coefficients,
    Eigen::Index draw,
    Eigen::Ref<Eigen::VectorXd> target) -> void
{
    assert(coefficients.rows() == projection.cols());
    assert(target.size() == projection.rows());
    target.setZero();
    for (CscReader::sparse_map_type<double>::InnerIterator it(
             coefficients, draw);
         it;
         ++it)
    {
        projection.axpy(
            static_cast<Eigen::Index>(it.row()), it.value(), target);
    }
}

}  // namespace gelex
