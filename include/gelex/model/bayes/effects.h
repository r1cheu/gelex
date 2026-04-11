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

#ifndef GELEX_MODEL_BAYES_EFFECTS_H_
#define GELEX_MODEL_BAYES_EFFECTS_H_

#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

struct RandomEffect
{
    RandomEffect(
        std::string name,
        std::vector<std::string> levels,
        Eigen::MatrixXd&& X)
        : name(std::move(name)), levels(std::move(levels)), X(std::move(X))
    {
        cols_squared_norm = this->X.colwise().squaredNorm();
    }

    std::string name;
    std::vector<std::string> levels;

    Eigen::MatrixXd X;
    Eigen::VectorXd cols_squared_norm;
};

struct GeneticEffect
{
    GeneticEffect(GeneticMode type, GenotypeStorage&& X)
        : type(type), X(std::move(X))
    {
        cols_squared_norm = get_matrix_ref(this->X).colwise().squaredNorm();
    }

    GeneticMode type;
    GenotypeStorage X;
    Eigen::VectorXd cols_squared_norm;

    auto is_monomorphic(Eigen::Index snp_index) const -> bool
    {
        return is_monomorphic_variant(X, snp_index);
    }

    auto num_mono() const -> Eigen::Index { return num_mono_variant(X); }
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_EFFECTS_H_
