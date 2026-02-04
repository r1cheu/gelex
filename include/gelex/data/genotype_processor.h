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

#ifndef GELEX_DATA_GENOTYPE_PROCESSOR_H_
#define GELEX_DATA_GENOTYPE_PROCESSOR_H_

#include <concepts>

#include <omp.h>
#include <Eigen/Dense>

namespace gelex
{

constexpr double GENOTYPE_EPSILON = 1e-10;

struct VariantStats
{
    double mean{0.0};
    double stddev{0.0};
    bool is_monomorphic{false};
};

template <typename T>
concept GenotypeProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process_variant(variant) } -> std::same_as<VariantStats>;
          { T::dom } -> std::convertible_to<bool>;
      };

// =========================================================================
// Method traits: each defines process_additive + process_dominant
// =========================================================================

struct CenterMethod
{
    static auto process_additive(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
    static auto process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
};

struct StandardizeMethod
{
    static auto process_additive(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
    static auto process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
};

struct OrthCenterMethod
{
    static auto process_additive(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
    static auto process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
};

struct OrthStandardizeMethod
{
    static auto process_additive(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
    static auto process_dominant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats;
};

// =========================================================================
// Template adapters: adapt Method traits to GenotypeProcessor concept
// =========================================================================

template <typename Method>
struct AdditiveProcessor
{
    static constexpr bool dom = false;
    static auto process_variant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        return Method::process_additive(variant);
    }
};

template <typename Method>
struct DominantProcessor
{
    static constexpr bool dom = true;
    static auto process_variant(Eigen::Ref<Eigen::VectorXd> variant)
        -> VariantStats
    {
        return Method::process_dominant(variant);
    }
};

// =========================================================================
// Batch processing template
// =========================================================================

template <GenotypeProcessor P>
auto process_matrix(
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs = nullptr) -> void
{
#pragma omp parallel for default(none) shared(genotype, freqs)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        auto col = genotype.col(i);
        auto stats = P::process_variant(col);
        if (freqs != nullptr)
        {
            (*freqs)(i) = stats.mean / 2;
        }
    }
}

// =========================================================================
// GRM method bundles
// =========================================================================

namespace grm
{

struct OrthStandardized
{
    using Additive = AdditiveProcessor<OrthStandardizeMethod>;
    using Dominant = DominantProcessor<OrthStandardizeMethod>;
};

struct Centered
{
    using Additive = AdditiveProcessor<CenterMethod>;
    using Dominant = DominantProcessor<CenterMethod>;
};

struct OrthCentered
{
    using Additive = AdditiveProcessor<OrthCenterMethod>;
    using Dominant = DominantProcessor<OrthCenterMethod>;
};

}  // namespace grm

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_PROCESSOR_H_
