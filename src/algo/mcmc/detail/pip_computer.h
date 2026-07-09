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

#ifndef GELEX_SRC_ALGO_MCMC_DETAIL_PIP_COMPUTER_H_
#define GELEX_SRC_ALGO_MCMC_DETAIL_PIP_COMPUTER_H_

#include <Eigen/Core>
#include <utility>

namespace gelex::detail
{

class PipComputer
{
   public:
    auto single(const Eigen::Ref<const Eigen::MatrixXd>& probabilities) const
        -> Eigen::VectorXd;
    auto joint(const Eigen::Ref<const Eigen::MatrixXd>& probabilities) const
        -> std::pair<Eigen::VectorXd, Eigen::VectorXd>;
};

}  // namespace gelex::detail

#endif  // GELEX_SRC_ALGO_MCMC_DETAIL_PIP_COMPUTER_H_
