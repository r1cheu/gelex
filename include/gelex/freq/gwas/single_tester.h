// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_SINGLE_TESTER_H_
#define GELEX_FREQ_GWAS_SINGLE_TESTER_H_

#include <Eigen/Core>

#include "gelex/data/genotype_method.h"
#include "gelex/freq/gwas/assoc_output.h"
#include "gelex/freq/gwas/assoc_tester.h"
#include "gelex/freq/reml/operators.h"

namespace gelex
{

class SingleTester final : public AssocTester
{
   public:
    SingleTester(GeneticMode mode, GenotypeMethod method);

    auto resize(Eigen::Index n_samples, Eigen::Index chunk_size)
        -> void override;

    [[nodiscard]] auto run(
        const LocusEncoder& encoder,
        Eigen::Index start,
        const GwasOperators& reml) -> TestResults override;

   private:
    static auto wald_test(
        Eigen::Ref<Eigen::MatrixXd> Z,
        Eigen::Ref<Eigen::MatrixXd> W,
        const GwasOperators& reml,
        AssocOutput& output) -> void;

    GeneticMode mode_;
    GenotypeMethod method_;

    Eigen::MatrixXd Z_;
    Eigen::MatrixXd W_;
    Eigen::VectorXd freqs_;
    AssocOutput output_;
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_SINGLE_TESTER_H_
