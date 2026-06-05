/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_ALGO_GWAS_ASSOC_TESTER_H_
#define GELEX_ALGO_GWAS_ASSOC_TESTER_H_

#include <memory>
#include <optional>
#include <span>

#include <Eigen/Core>

#include "gelex/algo/gwas/assoc_type.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

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

    [[nodiscard]] virtual auto genotype_buffer()
        -> Eigen::Ref<Eigen::MatrixXd> = 0;

    [[nodiscard]] virtual auto run(const RemlResult& reml) -> TestResults = 0;

    [[nodiscard]] static auto make(
        AssocType type,
        GeneticMode mode,
        GenotypeProcessMethod geno_method) -> std::unique_ptr<AssocTester>;

   protected:
    AssocTester() = default;
};
}  // namespace gelex

#endif  // GELEX_ALGO_GWAS_ASSOC_TESTER_H_
