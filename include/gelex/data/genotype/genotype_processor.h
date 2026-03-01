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

#include <type_traits>

#include <fmt/base.h>
#include <omp.h>
#include <Eigen/Dense>

#include "gelex/exception.h"
#include "gelex/internal/genotype_processor/genotype_processor.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{
template <typename T>
concept GenotypeProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process(variant) } -> std::same_as<LocusStatistic>;
      };

enum class GenotypeProcessMethod : uint8_t
{
    StandardizeHWE = 1,
    CenterHWE,
    OrthStandardizeHWE,
    OrthCenterHWE,
    Standardize,
    Center,
    OrthStandardize,
    OrthCenter
};

template <GeneticEffectType GT, detail::StatisticPolicy Stats, bool Scale>
using RawProcessor
    = detail::GenotypeProcessorStrategy<detail::RawPolicy<GT>, Stats, Scale>;

template <GeneticEffectType GT>
using Standardize = RawProcessor<GT, detail::SamplePolicy<GT>, true>;
template <GeneticEffectType GT>
using Center = RawProcessor<GT, detail::SamplePolicy<GT>, false>;

template <GeneticEffectType GT>
using StandardizeHWE = RawProcessor<GT, detail::HWEPolicy<GT>, true>;
template <GeneticEffectType GT>
using CenterHWE = RawProcessor<GT, detail::HWEPolicy<GT>, false>;

template <GeneticEffectType GT, detail::StatisticPolicy Stats, bool Scale>
using OrthProcessor = detail::
    GenotypeProcessorStrategy<detail::OrthogonalPolicy<GT>, Stats, Scale>;

template <GeneticEffectType GT>
using OrthStandardize = OrthProcessor<GT, detail::SamplePolicy<GT>, true>;
template <GeneticEffectType GT>
using OrthCenter = OrthProcessor<GT, detail::SamplePolicy<GT>, false>;

template <GeneticEffectType GT>
using OrthStandardizeHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, true>;
template <GeneticEffectType GT>
using OrthCenterHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, false>;

template <typename T, GeneticEffectType GT>
constexpr bool is_center_method_v
    = std::is_same_v<T, Center<GT>> || std::is_same_v<T, CenterHWE<GT>>
      || std::is_same_v<T, OrthCenter<GT>>
      || std::is_same_v<T, OrthCenterHWE<GT>>;

inline auto is_center_family_method(GenotypeProcessMethod method) -> bool
{
    switch (method)
    {
        case GenotypeProcessMethod::CenterHWE:
        case GenotypeProcessMethod::OrthCenterHWE:
        case GenotypeProcessMethod::Center:
        case GenotypeProcessMethod::OrthCenter:
            return true;
        case GenotypeProcessMethod::StandardizeHWE:
        case GenotypeProcessMethod::OrthStandardizeHWE:
        case GenotypeProcessMethod::Standardize:
        case GenotypeProcessMethod::OrthStandardize:
            return false;
    }
    throw InvalidInputException("Invalid genotype process method.");
}

template <GeneticEffectType GT>
auto get_genotype_process_method(GenotypeProcessMethod method)
    -> LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
{
    switch (method)
    {
        case GenotypeProcessMethod::StandardizeHWE:
            return &StandardizeHWE<GT>::process;
        case GenotypeProcessMethod::CenterHWE:
            return &CenterHWE<GT>::process;
        case GenotypeProcessMethod::OrthStandardizeHWE:
            return &OrthStandardizeHWE<GT>::process;
        case GenotypeProcessMethod::OrthCenterHWE:
            return &OrthCenterHWE<GT>::process;
        case GenotypeProcessMethod::Standardize:
            return &Standardize<GT>::process;
        case GenotypeProcessMethod::Center:
            return &Center<GT>::process;
        case GenotypeProcessMethod::OrthStandardize:
            return &OrthStandardize<GT>::process;
        case GenotypeProcessMethod::OrthCenter:
            return &OrthCenter<GT>::process;
    }
    throw InvalidInputException("Invalid genotype process method.");
}

template <GeneticEffectType GT>
auto get_center_genotype_method(GenotypeProcessMethod method)
    -> LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
{
    if (!is_center_family_method(method))
    {
        throw InvalidInputException(
            "assoc --geno-method supports only center-family methods: "
            "2 (center-hwe), 4 (orth-center-hwe), 6 (center), 8 (orth-center)");
    }
    return get_genotype_process_method<GT>(method);
}

template <GeneticEffectType GT>
auto process_matrix(
    GenotypeProcessMethod method,
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs = nullptr) -> void
{
    auto process = get_genotype_process_method<GT>(method);
#pragma omp parallel for default(none) shared(genotype, freqs, process)
    for (Eigen::Index i = 0; i < genotype.cols(); ++i)
    {
        auto col = genotype.col(i);
        auto stats = process(col);
        if (freqs != nullptr)
        {
            (*freqs)(i) = stats.maf;
        }
    }
}
}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::GenotypeProcessMethod> : formatter<string_view>
{
    auto format(gelex::GenotypeProcessMethod t, format_context& ctx) const
        -> format_context::iterator
    {
        string_view name = "unknown";
        using gelex::GenotypeProcessMethod;
        switch (t)
        {
            case GenotypeProcessMethod::StandardizeHWE:
                name = "StandardizeHWE";
                break;
            case GenotypeProcessMethod::CenterHWE:
                name = "CenterHWE";
                break;
            case GenotypeProcessMethod::OrthStandardizeHWE:
                name = "OrthStandardizeHWE";
                break;
            case GenotypeProcessMethod::OrthCenterHWE:
                name = "OrthCenterHWE";
                break;
            case GenotypeProcessMethod::Standardize:
                name = "Standardize";
                break;
            case GenotypeProcessMethod::Center:
                name = "Center";
                break;
            case GenotypeProcessMethod::OrthStandardize:
                name = "OrthStandardize";
                break;
            case GenotypeProcessMethod::OrthCenter:
                name = "OrthCenter";
                break;
            default:
                name = "Unknown Genotype Process Method";
                break;
        }

        return formatter<string_view>::format(name, ctx);
    }
};
}  // namespace fmt

#endif  // GELEX_DATA_GENOTYPE_PROCESSOR_H_
