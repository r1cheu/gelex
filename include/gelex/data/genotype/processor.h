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

#include "gelex/data/genotype/detail/processor_strategy.h"
#include "gelex/data/genotype/process_method.h"
#include "gelex/exception.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex::genotype
{
template <typename T>
concept GenotypeProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process(variant) } -> std::same_as<gelex::LocusStatistic>;
      };

template <gelex::GeneticMode GT, detail::StatisticPolicy Stats, bool Scale>
using RawProcessor
    = detail::ProcessorStrategy<detail::RawPolicy<GT>, Stats, Scale>;

template <gelex::GeneticMode GT>
using Standardize = RawProcessor<GT, detail::SamplePolicy<GT>, true>;
template <gelex::GeneticMode GT>
using Center = RawProcessor<GT, detail::SamplePolicy<GT>, false>;

template <gelex::GeneticMode GT>
using StandardizeHWE = RawProcessor<GT, detail::HWEPolicy<GT>, true>;
template <gelex::GeneticMode GT>
using CenterHWE = RawProcessor<GT, detail::HWEPolicy<GT>, false>;

template <gelex::GeneticMode GT, detail::StatisticPolicy Stats, bool Scale>
using OrthProcessor
    = detail::ProcessorStrategy<detail::OrthogonalPolicy<GT>, Stats, Scale>;

template <gelex::GeneticMode GT>
using OrthStandardize = OrthProcessor<GT, detail::SamplePolicy<GT>, true>;
template <gelex::GeneticMode GT>
using OrthCenter = OrthProcessor<GT, detail::SamplePolicy<GT>, false>;

template <gelex::GeneticMode GT>
using OrthStandardizeHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, true>;
template <gelex::GeneticMode GT>
using OrthCenterHWE = OrthProcessor<GT, detail::OrthHWEPolicy<GT>, false>;

template <gelex::GeneticMode GT, bool Scale>
using NOIAProcessor = detail::ProcessorStrategy<
    detail::IdentityPolicy<GT>,
    detail::NOIAPolicy<GT>,
    Scale>;

template <gelex::GeneticMode GT>
using NOIAStandardize = NOIAProcessor<GT, true>;
template <gelex::GeneticMode GT>
using NOIACenter = NOIAProcessor<GT, false>;

template <typename T, gelex::GeneticMode GT>
constexpr bool is_center_method_v
    = std::is_same_v<T, Center<GT>> || std::is_same_v<T, CenterHWE<GT>>
      || std::is_same_v<T, OrthCenter<GT>>
      || std::is_same_v<T, OrthCenterHWE<GT>>
      || std::is_same_v<T, NOIACenter<GT>>;

template <gelex::GeneticMode GT>
auto get_genotype_process_method(gelex::GenotypeProcessMethod method)
    -> gelex::LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
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

template <gelex::GeneticMode GT>
auto get_center_genotype_method(gelex::GenotypeProcessMethod method)
    -> gelex::LocusStatistic (*)(Eigen::Ref<Eigen::VectorXd>)
{
    if (!method.is_center())
    {
        throw gelex::GelexException(
            "assoc --geno-method supports only center-family methods: "
            "CenterHWE(CH), OrthCenterHWE(OCH), Center(C), OrthCenter(OC), "
            "NOIACenter(NC)");
    }
    return get_genotype_process_method<GT>(method);
}

template <gelex::GeneticMode GT>
auto process_matrix(
    gelex::GenotypeProcessMethod method,
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
}  // namespace gelex::genotype

#endif  // GELEX_DATA_GENOTYPE_PROCESSOR_H_
