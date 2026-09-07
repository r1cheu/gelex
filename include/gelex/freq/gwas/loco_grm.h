// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_LOCO_GRM_H_
#define GELEX_FREQ_GWAS_LOCO_GRM_H_

#include <Eigen/Core>

namespace gelex
{

class LocoGrmBuilder
{
   public:
    explicit LocoGrmBuilder(const Eigen::MatrixXd& whole_grm);

    auto build_into(
        const Eigen::Ref<const Eigen::MatrixXd>& chromosome_grm,
        Eigen::MatrixXd& target) const -> void;

   private:
    // Borrowed; the whole-genome GRM must outlive this builder.
    const Eigen::MatrixXd* whole_grm_;
    double whole_denominator_ = 0.0;
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_LOCO_GRM_H_
