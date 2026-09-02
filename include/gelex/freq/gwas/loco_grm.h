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
