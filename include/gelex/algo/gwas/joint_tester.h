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

#ifndef GELEX_ALGO_GWAS_JOINT_TESTER_H_
#define GELEX_ALGO_GWAS_JOINT_TESTER_H_

#include <Eigen/Core>

#include "gelex/algo/gwas/assoc_output.h"
#include "gelex/algo/gwas/assoc_tester.h"
#include "gelex/algo/reml/result.h"
#include "gelex/data/genotype/process_method.h"

namespace gelex
{

class JointTester final : public AssocTester
{
   public:
    explicit JointTester(GenotypeProcessMethod method);

    auto resize(Eigen::Index n_samples, Eigen::Index chunk_size)
        -> void override;

    [[nodiscard]] auto genotype_buffer()
        -> Eigen::Ref<Eigen::MatrixXd> override;

    [[nodiscard]] auto run(const RemlResult& reml) -> TestResults override;

   private:
    GenotypeProcessMethod method_;

    Eigen::MatrixXd raw_;
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

#endif  // GELEX_ALGO_GWAS_JOINT_TESTER_H_
