// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_JOINT_TESTER_H_
#define GELEX_FREQ_GWAS_JOINT_TESTER_H_

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/freq/gwas/assoc_output.h"
#include "gelex/freq/gwas/assoc_tester.h"
#include "gelex/freq/reml/operators.h"

namespace gelex
{

class JointTester final : public AssocTester
{
   public:
    explicit JointTester(GenotypeMethod method);

    auto resize(Eigen::Index n_samples, Eigen::Index chunk_size)
        -> void override;

    [[nodiscard]] auto run(
        const LocusEncoder& encoder,
        Eigen::Index start,
        const GwasOperators& reml) -> TestResults override;

   private:
    GenotypeMethod method_;

    Eigen::MatrixXd Z_a_;
    Eigen::MatrixXd Z_d_;
    Eigen::MatrixXd W_;
    Eigen::VectorXd freqs_;

    AssocOutput add_;
    AssocOutput dom_;
    Eigen::VectorXd zt_a_Pzd_;
    Eigen::VectorXd joint_p_;
    Eigen::VectorXd total_pve_;
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_JOINT_TESTER_H_
