// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_FREQ_GWAS_ASSOC_TESTER_H_
#define GELEX_FREQ_GWAS_ASSOC_TESTER_H_

#include <Eigen/Core>
#include <memory>
#include <optional>
#include <span>

#include "gelex/data/genotype_method.h"
#include "gelex/freq/gwas/assoc_type.h"
#include "gelex/freq/reml/operators.h"
#include "gelex/genetic_mode.h"

namespace gelex
{

class LocusEncoder;

struct TestResult
{
    std::span<const double> beta;
    std::span<const double> se;
    std::span<const double> p;
    std::span<const double> pve;
};

struct TestResults
{
    std::span<const double> freq;
    TestResult additive;
    std::optional<TestResult> dominance;
    std::optional<std::span<const double>> joint_p;
    std::optional<std::span<const double>> total_pve;
};

class AssocTester
{
   public:
    virtual ~AssocTester() = default;

    AssocTester(const AssocTester&) = delete;
    auto operator=(const AssocTester&) -> AssocTester& = delete;
    AssocTester(AssocTester&&) = delete;
    auto operator=(AssocTester&&) -> AssocTester& = delete;

    virtual auto resize(Eigen::Index n_samples, Eigen::Index chunk_size) -> void
        = 0;

    // Fuses the chunk [start, start + chunk_size) straight from its packed form
    // via encoder, then tests. chunk_size is the width fixed by the last
    // resize.
    [[nodiscard]] virtual auto run(
        const LocusEncoder& encoder,
        Eigen::Index start,
        const GwasOperators& reml) -> TestResults = 0;

    [[nodiscard]] static auto make(
        AssocType type,
        GeneticMode mode,
        GenotypeMethod geno_method) -> std::unique_ptr<AssocTester>;

   protected:
    AssocTester() = default;
};
}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_ASSOC_TESTER_H_
