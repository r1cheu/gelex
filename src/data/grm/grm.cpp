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

#include "gelex/data/grm/grm.h"

#include <algorithm>
#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <Eigen/Core>
#ifdef USE_MKL
#include <mkl.h>
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include "gelex/data/bed.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/infra/logging/notify.h"

namespace gelex
{
using Eigen::Index;

namespace
{

auto update_grm(
    Eigen::Ref<Eigen::MatrixXd> grm,
    const Eigen::Ref<const Eigen::MatrixXd>& genotype) -> void
{
    const auto n = static_cast<int>(genotype.rows());
    const auto m = static_cast<int>(genotype.cols());

    // dsyrk: C := alpha * A * A^T + beta * C, CblasLower fills lower triangle.
    cblas_dsyrk(
        CblasColMajor,
        CblasLower,
        CblasNoTrans,
        n,
        m,
        1.0,
        genotype.data(),
        n,
        1.0,
        grm.data(),
        n);
}

}  // namespace

auto chromosome_ranges(const dataframe::DataFrame<std::string>& bim)
    -> std::vector<GrmRange>
{
    const auto num_snps = static_cast<Index>(bim.rows());
    auto chrom = bim["chrom"].as<std::string>();

    std::vector<GrmRange> ranges;
    std::string current;
    Index range_start = 0;
    for (Index i = 0; i < num_snps; ++i)
    {
        if (chrom[static_cast<std::size_t>(i)] != current)
        {
            if (!current.empty())
            {
                ranges.push_back({current, range_start, i});
            }
            current = chrom[static_cast<std::size_t>(i)];
            range_start = i;
        }
    }
    if (!current.empty())
    {
        ranges.push_back({current, range_start, num_snps});
    }
    return ranges;
}

GrmBuilder::GrmBuilder(
    const Bed& bed,
    GeneticModeSet modes,
    GenotypeMethod method,
    Index chunk_size,
    GrmObserver observer)
    : bed_(bed),
      modes_(modes),
      method_(method),
      chunk_size_(chunk_size),
      observer_(std::move(observer))
{
}

auto GrmBuilder::accumulate(std::string_view label, Index start, Index end)
    -> std::vector<GrmMatrix>
{
    const Index n = bed_.num_samples();
    std::vector<Eigen::MatrixXd> grms(
        modes_.size(), Eigen::MatrixXd::Zero(n, n));

    for (Index s = start; s < end; s += chunk_size_)
    {
        const Index e = std::min(s + chunk_size_, end);
        const Eigen::MatrixXd genotype = bed_.read<double>(s, e);
        Eigen::MatrixXd encoded(n, e - s);

        std::size_t k = 0;
        for (GeneticMode mode : modes_.each())
        {
            encode_into<double, double>(genotype, encoded, mode, method_);
            update_grm(grms[k], encoded);
            ++k;
        }

        processed_ += (e - s);
        notify(
            observer_,
            GrmProgressEvent{
                static_cast<size_t>(processed_),
                static_cast<size_t>(total_),
                false});
    }

    std::vector<GrmMatrix> results;
    results.reserve(grms.size());
    std::size_t k = 0;
    for (GeneticMode mode : modes_.each())
    {
        const double denominator = grms[k].trace() / static_cast<double>(n);
        results.push_back(
            {std::string{label}, mode, std::move(grms[k]), denominator});
        ++k;
    }
    return results;
}

auto GrmBuilder::build(std::span<const GrmRange> ranges, const Sink& sink)
    -> void
{
    total_ = 0;
    for (const auto& range : ranges)
    {
        total_ += (range.end - range.start);
    }
    processed_ = 0;

    for (const auto& range : ranges)
    {
        for (const auto& matrix :
             accumulate(range.label, range.start, range.end))
        {
            sink(matrix);
        }
    }

    notify(
        observer_,
        GrmProgressEvent{
            static_cast<size_t>(processed_),
            static_cast<size_t>(total_),
            true});
}

}  // namespace gelex
