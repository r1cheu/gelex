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

#include "algo/mcmc/detail/pip_computer.h"

#include <Eigen/Core>
#include <utility>

namespace gelex::detail
{

auto PipComputer::single(const Eigen::Ref<const Eigen::MatrixXd>& probabilities)
    const -> Eigen::VectorXd
{
    return (Eigen::VectorXd::Ones(probabilities.rows()) - probabilities.col(0))
        .eval();
}

auto PipComputer::joint(const Eigen::Ref<const Eigen::MatrixXd>& probabilities)
    const -> std::pair<Eigen::VectorXd, Eigen::VectorXd>
{
    auto additive = (probabilities.col(1) + probabilities.col(3)).eval();
    auto dominance = (probabilities.col(2) + probabilities.col(3)).eval();
    return {std::move(additive), std::move(dominance)};
}

}  // namespace gelex::detail
