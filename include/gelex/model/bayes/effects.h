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

#include "gelex/data/genotype/genotype.h"
#include "gelex/infra/stats/detail/var.h"
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
        XtX_diag = this->X.colwise().squaredNorm();
    }

    std::string name;
    std::vector<std::string> levels;

    Eigen::MatrixXd X;
    Eigen::VectorXd XtX_diag;
};

struct GeneticEffect
{
    GeneticEffect(GeneticMode type, gelex::genotype::Genotype&& X)
        : type(type), X(std::move(X))
    {
        XtX_diag = this->X.matrix().colwise().squaredNorm();
        variance = gelex::stats::detail::var<0>(X.matrix()).sum();
    }

    GeneticMode type;
    gelex::genotype::Genotype X;
    Eigen::VectorXd XtX_diag;
    double variance{};

    auto is_monomorphic(Eigen::Index snp_index) const -> bool
    {
        return X.is_monomorphic(snp_index);
    }

    auto num_mono() const -> Eigen::Index { return X.num_mono(); }
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_EFFECTS_H_
