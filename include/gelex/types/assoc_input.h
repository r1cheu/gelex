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

#ifndef GELEX_TYPES_ASSOC_INPUT_H_
#define GELEX_TYPES_ASSOC_INPUT_H_

#include <Eigen/Core>

namespace gelex
{

struct AssocInput
{
    AssocInput() = default;

    AssocInput(Eigen::Index n_samples, Eigen::Index chunk_size)
    {
        resize(n_samples, chunk_size);
    }

    auto resize(Eigen::Index n_samples, Eigen::Index chunk_size) -> void
    {
        Z.resize(n_samples, chunk_size);
        W.resize(n_samples, chunk_size);
        V_inv.resize(n_samples, n_samples);
        V_inv_y.resize(n_samples);
    }

    Eigen::MatrixXd Z;      // SNP matrix
    Eigen::MatrixXd V_inv;  // Inverse of covariance matrix
    Eigen::VectorXd V_inv_y;
    Eigen::MatrixXd W;  // Intermediate buffer for V^{-1} Z
};

struct AssocOutput
{
    AssocOutput() = default;

    explicit AssocOutput(Eigen::Index chunk_size) { resize(chunk_size); }

    auto resize(Eigen::Index chunk_size) -> void
    {
        beta.resize(chunk_size);
        se.resize(chunk_size);
        stats.resize(chunk_size);
        p_value.resize(chunk_size);
        zt_v_inv_r.resize(chunk_size);
        zt_v_inv_z.resize(chunk_size);
    }

    Eigen::VectorXd beta;
    Eigen::VectorXd se;
    Eigen::VectorXd stats;
    Eigen::VectorXd p_value;

    Eigen::VectorXd zt_v_inv_r;  // Z^t V^{-1} residual
    Eigen::VectorXd zt_v_inv_z;  // Z^t V^{-1} Z
};

}  // namespace gelex

#endif  // GELEX_TYPES_ASSOC_INPUT_H_
