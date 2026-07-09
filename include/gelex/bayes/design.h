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

#ifndef GELEX_BAYES_DESIGN_H_
#define GELEX_BAYES_DESIGN_H_

#include <Eigen/Core>
#include <cstdint>
#include <string>
#include <vector>

#include "gelex/data/genotype.h"
#include "gelex/infra/stats/detail/var.h"
#include "gelex/types/genetic_mode.h"

namespace gelex
{
class FieldVisitor;
}

namespace gelex::bayes
{

struct RandomDesign
{
    RandomDesign(
        std::string name,
        std::vector<std::string> levels,
        Eigen::MatrixXd&& X)
        : name(std::move(name)), levels(std::move(levels)), X(std::move(X))
    {
        XtX_diag = this->X.colwise().squaredNorm();
    }

    std::string name;
    std::vector<std::string> levels;

    Eigen::MatrixXd X;
    Eigen::VectorXd XtX_diag;

    auto visit(FieldVisitor& visitor) const -> void;
};

struct GeneticDesign
{
    GeneticDesign(GeneticMode type, gelex::Genotype&& X)
        : type(type), X(std::move(X))
    {
        XtX_diag = this->X.matrix().colwise().squaredNorm();
        col_var = gelex::detail::matvar<0>(
            this->X.matrix(), gelex::detail::VarNormType::Population);
    }

    GeneticMode type;
    gelex::Genotype X;
    Eigen::VectorXd XtX_diag;
    Eigen::RowVectorXd col_var;

    auto valid_indices() const -> const std::vector<int64_t>&
    {
        return X.valid_indices();
    }

    auto num_invalid() const -> Eigen::Index { return X.num_invalid(); }
};

}  // namespace gelex::bayes

#endif  // GELEX_BAYES_DESIGN_H_
