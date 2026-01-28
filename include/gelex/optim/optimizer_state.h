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

#ifndef GELEX_OPTIM_OPTIMIZER_STATE_H_
#define GELEX_OPTIM_OPTIMIZER_STATE_H_

#include <Eigen/Core>

namespace gelex
{

class FreqModel;

class OptimizerState
{
   public:
    explicit OptimizerState(const FreqModel& model);

    auto phenotype_variance() const -> double { return phenotype_variance_; }
    auto num_individuals() const -> Eigen::Index { return num_individuals_; }

    // computed matrices, public for Policy access
    Eigen::MatrixXd v;
    Eigen::MatrixXd proj;
    Eigen::VectorXd proj_y;
    Eigen::MatrixXd tx_vinv_x;
    double logdet_v{};

    // for AI policy
    Eigen::MatrixXd hess_inv;
    Eigen::MatrixXd dvpy;        // n x n_comp, each column is K_i * Py
    Eigen::VectorXd first_grad;  // first derivative

   private:
    Eigen::Index num_individuals_{};
    double phenotype_variance_{};
};

}  // namespace gelex

#endif  // GELEX_OPTIM_OPTIMIZER_STATE_H_
