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

#ifndef GELEX_FREQ_GWAS_ASSOC_OUTPUT_H_
#define GELEX_FREQ_GWAS_ASSOC_OUTPUT_H_

#include <Eigen/Core>

namespace gelex
{

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
        pve.resize(chunk_size);
        zt_Pr.resize(chunk_size);
        zt_Pz.resize(chunk_size);
    }

    Eigen::VectorXd beta;
    Eigen::VectorXd se;
    Eigen::VectorXd stats;
    Eigen::VectorXd p_value;
    Eigen::VectorXd pve;

    Eigen::VectorXd zt_Pr;
    Eigen::VectorXd zt_Pz;
};

}  // namespace gelex

#endif  // GELEX_FREQ_GWAS_ASSOC_OUTPUT_H_
