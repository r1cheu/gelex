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

#ifndef GELEX_SRC_ALGO_MCMC_DETAIL_PVE_COMPUTER_H_
#define GELEX_SRC_ALGO_MCMC_DETAIL_PVE_COMPUTER_H_

#include <Eigen/Core>

#include "gelex/types/genetic_mode.h"

namespace gelex
{

namespace bayes
{
class GeneticDesign;
}

namespace detail
{

class PveComputer
{
   public:
    PveComputer(const bayes::GeneticDesign& design, double phenotype_var);

    auto single(GeneticMode mode, const Eigen::Ref<const Eigen::VectorXd>& beta)
        const -> Eigen::VectorXd;
    auto total(
        const Eigen::Ref<const Eigen::VectorXd>& beta_A,
        const Eigen::Ref<const Eigen::VectorXd>& beta_D) const
        -> Eigen::VectorXd;

   private:
    const bayes::GeneticDesign& design_;
    double phenotype_var_{};
    Eigen::RowVectorXd cov_ad_;
};

}  // namespace detail

}  // namespace gelex

#endif  // GELEX_SRC_ALGO_MCMC_DETAIL_PVE_COMPUTER_H_
