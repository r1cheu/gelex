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

#include <optional>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/infra/stats/conjugate_prior.h"
#include "gelex/model/bayes/genotype_storage.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

struct MixtureConfig
{
    Eigen::VectorXd init_proportion;
    std::optional<Eigen::VectorXd> scale;
    bool estimate_pi{false};
};

struct RandomEffect
{
    RandomEffect(
        std::optional<std::vector<std::string>> levels,
        Eigen::MatrixXd&& X)
        : X(std::move(X)), levels(std::move(levels))
    {
        cols_squared_norm = this->X.colwise().squaredNorm();
    }

    Eigen::MatrixXd X;
    Eigen::VectorXd cols_squared_norm;

    std::optional<std::vector<std::string>> levels;

    detail::ScaledInvChiSqParams prior{4, 0};
    double init_variance{0.0};
};

struct SignConfig
{
    double init_positive_prob{0.5};
};

struct GeneticEffect
{
    GeneticEffect(GeneticEffectType type, GenotypeStorage&& X)
        : type(type), X(std::move(X))
    {
        cols_squared_norm = get_matrix_ref(this->X).colwise().squaredNorm();
    }

    GeneticEffectType type;
    GenotypeStorage X;
    Eigen::VectorXd cols_squared_norm;

    detail::ScaledInvChiSqParams marker_variance_prior{4, 0};
    double init_marker_variance{0.0};
    Eigen::Index marker_variance_size{0};

    std::optional<MixtureConfig> mixture;
    std::optional<SignConfig> sign;

    auto is_monomorphic(Eigen::Index snp_index) const -> bool
    {
        return is_monomorphic_variant(X, snp_index);
    }

    auto num_mono() const -> Eigen::Index { return num_mono_variant(X); }
};

struct Residual
{
    detail::ScaledInvChiSqParams prior{-2, 0};
    double init_variance{0.0};
};

}  // namespace gelex::bayes

#endif  // GELEX_MODEL_BAYES_EFFECTS_H_
