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
#ifndef GELEX_UTILS_PHENOTYPE_TRANSFORMER_H
#define GELEX_UTILS_PHENOTYPE_TRANSFORMER_H

#include <Eigen/Dense>

namespace gelex::detail
{

class PhenotypeTransformer
{
   public:
    explicit PhenotypeTransformer(double offset = 3.0 / 8.0);

    auto apply_dint(Eigen::Ref<Eigen::VectorXd> phenotype) const -> void;
    auto apply_iint(
        Eigen::Ref<Eigen::VectorXd> phenotype,
        const Eigen::Ref<const Eigen::MatrixXd>& covariates) const -> void;

   private:
    double offset_;

    auto int_transform(Eigen::Ref<Eigen::VectorXd> values) const -> void;
    static auto compute_ranks(const Eigen::Ref<const Eigen::VectorXd>& values)
        -> Eigen::VectorXd;
    static auto compute_residuals(
        const Eigen::Ref<const Eigen::VectorXd>& y,
        const Eigen::Ref<const Eigen::MatrixXd>& X) -> Eigen::VectorXd;
};

}  // namespace gelex::detail

#endif  // GELEX_UTILS_PHENOTYPE_TRANSFORMER_H
