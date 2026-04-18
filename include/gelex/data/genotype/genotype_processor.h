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
#include "gelex/types/genotype_process_method.h"

namespace gelex
{
template <typename T>
concept GenotypeProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process(variant) } -> std::same_as<LocusStatistic>;
      };

template <GeneticMode GT, detail::StatisticPolicy Stats, bool Scale>
using RawProcessor
    = detail::GenotypeProcessorStrategy<detail::RawPolicy<GT>, Stats, Scale>;

template <GeneticMode GT>
using Standardize = RawProcessor<GT, detail::SamplePolicy<GT>, true>;
template <GeneticMode GT>
using Center = RawProcessor<GT, detail::SamplePolicy<GT>, false>;

template <GeneticMode GT>
using StandardizeHWE = RawProcessor<GT, detail::HWEPolicy<GT>, true>;
template <GeneticMode GT>
using CenterHWE = RawProcessor<GT, detail::HWEPolicy<GT>, false>;

template <GeneticMode GT, detail::StatisticPolicy Stats, bool Scale>
using OrthProcessor = detail::
    GenotypeProcessorStrategy<detail::OrthogonalPolicy<GT>, Stats, Scale>;

template <GeneticMode GT>
using OrthStandardize = OrthProcessor<GT, detail::SamplePolicy<GT>, true>;
template <GeneticMode GT>
using OrthCenter = OrthProcessor<GT, detail::SamplePolicy<GT>, false>;

template <GeneticMode GT>
using OrthStandardizeHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, true>;
template <GeneticMode GT>
using OrthCenterHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, false>;

template <GeneticMode GT, bool Scale>
using NOIAProcessor = detail::GenotypeProcessorStrategy<
    detail::IdentityPolicy<GT>,
    detail::NOIAPolicy<GT>,
    Scale>;

template <GeneticMode GT>
using NOIAStandardize = NOIAProcessor<GT, true>;
template <GeneticMode GT>
using NOIACenter = NOIAProcessor<GT, false>;

template <typename T, GeneticMode GT>
constexpr bool is_center_method_v
    = std::is_same_v<T, Center<GT>> || std::is_same_v<T, CenterHWE<GT>>
      || std::is_same_v<T, OrthCenter<GT>>
      || std::is_same_v<T, OrthCenterHWE<GT>>
      || std::is_same_v<T, NOIACenter<GT>>;

template <GeneticMode GT>
auto get_genotype_process_method(GenotypeProcessMethod method)
    -> LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
{
    if (method.is_noia())
    {
        return method.is_center() ? &NOIACenter<GT>::process
                                  : &NOIAStandardize<GT>::process;
    }
    if (method.is_orthogonal())
    {
        if (method.is_hwe())
        {
            return method.is_center() ? &OrthCenterHWE<GT>::process
                                      : &OrthStandardizeHWE<GT>::process;
        }
        return method.is_center() ? &OrthCenter<GT>::process
                                  : &OrthStandardize<GT>::process;
    }
    if (method.is_hwe())
    {
        return method.is_center() ? &CenterHWE<GT>::process
                                  : &StandardizeHWE<GT>::process;
    }
    return method.is_center() ? &Center<GT>::process
                              : &Standardize<GT>::process;
}

template <GeneticMode GT>
auto get_center_genotype_method(GenotypeProcessMethod method)
    -> LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
{
    if (!method.is_center())
    {
        throw GelexException(
            "assoc --geno-method supports only center-family methods: "
            "CenterHWE(CH), OrthCenterHWE(OCH), Center(C), OrthCenter(OC), "
            "NOIACenter(NC)");
    }
    return get_genotype_process_method<GT>(method);
}

template <GeneticMode GT>
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

#endif  // GELEX_DATA_GENOTYPE_PROCESSOR_H_
