// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/freq/gwas/loco_grm.h"

#include <Eigen/Core>
#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex
{

LocoGrmBuilder::LocoGrmBuilder(const Eigen::MatrixXd& whole_grm)
    : whole_grm_(&whole_grm)
{
    if (whole_grm.rows() == 0 || whole_grm.rows() != whole_grm.cols())
    {
        throw GelexException(
            fmt::format(
                "LOCO error: whole GRM must be a non-empty square matrix, got "
                "{}x{}",
                whole_grm.rows(),
                whole_grm.cols()));
    }
    whole_denominator_
        = whole_grm.trace() / static_cast<double>(whole_grm.rows());
}

auto LocoGrmBuilder::build_into(
    const Eigen::Ref<const Eigen::MatrixXd>& chromosome_grm,
    Eigen::MatrixXd& target) const -> void
{
    if (chromosome_grm.rows() != whole_grm_->rows()
        || chromosome_grm.cols() != whole_grm_->cols())
    {
        throw GelexException(
            fmt::format(
                "LOCO error: chromosome GRM shape {}x{} does not match whole "
                "GRM shape {}x{}",
                chromosome_grm.rows(),
                chromosome_grm.cols(),
                whole_grm_->rows(),
                whole_grm_->cols()));
    }

    const double chromosome_denominator
        = chromosome_grm.trace() / static_cast<double>(chromosome_grm.rows());
    const double loco_denominator = whole_denominator_ - chromosome_denominator;
    if (loco_denominator <= 0)
    {
        throw GelexException(
            fmt::format(
                "LOCO error: Chromosome GRM denominator ({}) is greater than "
                "or equal to Whole GRM denominator ({})",
                chromosome_denominator,
                whole_denominator_));
    }

    target = (*whole_grm_ - chromosome_grm) / loco_denominator;
}

}  // namespace gelex
