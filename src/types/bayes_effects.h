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

#ifndef GELEX_TYPES_BAYES_EFFECTS_H_
#define GELEX_TYPES_BAYES_EFFECTS_H_
#include <string>
#include <variant>
#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/distribution.h"
#include "../src/types/fixed_effects.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_mmap.h"

namespace gelex
{
namespace bayes
{

// Variant type for genotype storage
using GenotypeStorage = std::variant<GenotypeMap, GenotypeMatrix>;

// Helper to get matrix reference from either storage type
// Returns Eigen::Ref which can wrap both Map and MatrixXd
inline Eigen::Ref<const Eigen::MatrixXd> get_matrix_ref(
    const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> Eigen::Ref<const Eigen::MatrixXd>
        { return s.matrix(); },
        storage);
}

// Helper to get number of rows
inline Eigen::Index get_rows(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.rows(); }, storage);
}

// Helper to get number of columns
inline Eigen::Index get_cols(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.cols(); }, storage);
}

inline const Eigen::VectorXd& get_means(const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.mean(); },
        storage);
}

inline const Eigen::VectorXd& get_stddev(const GenotypeStorage& storage)
{
    return std::visit(
        [](const auto& s) -> const Eigen::VectorXd& { return s.stddev(); },
        storage);
}

// Helper to check monomorphic status
inline bool is_monomorphic_variant(
    const GenotypeStorage& storage,
    Eigen::Index idx)
{
    return std::visit(
        [idx](const auto& s) { return s.is_monomorphic(idx); }, storage);
}

// Helper to get number of monomorphic markers
inline Eigen::Index num_mono_variant(const GenotypeStorage& storage)
{
    return std::visit([](const auto& s) { return s.num_mono(); }, storage);
}

struct Pi
{
    Eigen::VectorXd prop;
    Eigen::VectorXi count;
};

struct FixedState
{
    explicit FixedState(const FixedEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())) {};
    Eigen::VectorXd coeffs;
};

struct RandomEffect
{
    RandomEffect(
        std::optional<std::vector<std::string>> levels,
        Eigen::MatrixXd&& X)
        : X(std::move(X)), levels(std::move(levels))
    {
        cols_norm = this->X.colwise().squaredNorm();
    }

    Eigen::MatrixXd X;
    Eigen::VectorXd cols_norm;

    std::optional<std::vector<std::string>> levels;

    detail::ScaledInvChiSqParams prior{4, 0};
    double init_variance{0.0};
};

struct RandomState
{
    explicit RandomState(const RandomEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(effect.X.cols())),
          variance{effect.init_variance}
    {
    }

    Eigen::VectorXd coeffs;
    double variance{0.0};
};

struct GeneticEffect
{
    explicit GeneticEffect(GenotypeMap&& X) : X(std::move(X))
    {
        cols_norm = get_matrix_ref(this->X).colwise().squaredNorm();
    }

    explicit GeneticEffect(GenotypeMatrix&& X) : X(std::move(X))
    {
        cols_norm = get_matrix_ref(this->X).colwise().squaredNorm();
    }

    explicit GeneticEffect(GenotypeStorage&& X) : X(std::move(X))
    {
        cols_norm = get_matrix_ref(this->X).colwise().squaredNorm();
    }

    GenotypeStorage X;
    Eigen::VectorXd cols_norm;

    detail::ScaledInvChiSqParams marker_variance_prior{4, 0};
    double init_marker_variance{0.0};
    Eigen::Index marker_variance_size{0};

    std::optional<Eigen::VectorXd> init_pi;
    std::optional<Eigen::VectorXd> scale;
    bool estimate_pi{false};

    bool is_monomorphic(Eigen::Index snp_index) const
    {
        return is_monomorphic_variant(X, snp_index);
    }

    Eigen::Index num_mono() const { return num_mono_variant(X); }
};

struct AdditiveEffect : GeneticEffect
{
    using GeneticEffect::GeneticEffect;
};

struct GeneticState
{
    explicit GeneticState(const GeneticEffect& effect)
        : coeffs(Eigen::VectorXd::Zero(bayes::get_cols(effect.X))),
          u(Eigen::VectorXd::Zero(bayes::get_rows(effect.X))),
          marker_variance(
              Eigen::VectorXd::Constant(
                  effect.marker_variance_size,
                  effect.init_marker_variance))

    {
        if (effect.init_pi)
        {
            tracker = Eigen::VectorXi::Zero(bayes::get_cols(effect.X));
            pi
                = {effect.init_pi.value(),
                   Eigen::VectorXi::Zero(effect.init_pi->size())};

            if (const auto num_components = effect.init_pi->size();
                num_components > 2)
            {
                component_u.resize(num_components - 1);
                for (auto& vec : component_u)
                {
                    vec = Eigen::VectorXd::Zero(u.size());
                }
                component_variance = Eigen::VectorXd::Zero(num_components - 1);
            }
        };
    }
    Eigen::VectorXd coeffs;
    Eigen::VectorXd u;

    Eigen::VectorXi tracker;
    Pi pi;

    double variance{};
    double heritability{};
    Eigen::VectorXd marker_variance;

    std::vector<Eigen::VectorXd> component_u;
    Eigen::VectorXd component_variance;
};

struct AdditiveState : GeneticState
{
    using GeneticState::GeneticState;
};

struct DominantEffect : GeneticEffect
{
    using GeneticEffect::GeneticEffect;
};

struct DominantState : GeneticState
{
    using GeneticState::GeneticState;
};

struct Residual
{
    detail::ScaledInvChiSqParams prior{-2, 0};
    double init_variance{0.0};
};

struct ResidualState
{
    Eigen::VectorXd y_adj;
    double variance{0.0};
};
}  // namespace bayes
}  // namespace gelex

#endif  // GELEX_TYPES_BAYES_EFFECTS_H_
