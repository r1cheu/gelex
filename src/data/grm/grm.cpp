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

#include <Eigen/Core>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gelex/data/bed.h"
#include "gelex/data/encode/encoder.h"
#include "gelex/data/encode/spec.h"
#include "gelex/data/encode/stats.h"
#include "gelex/data/encode/types.h"
#include "gelex/data/marker_range.h"
#include "gelex/infra/notify.h"

namespace gelex
{
using Eigen::Index;

GrmBuilder::GrmBuilder(
    const Bed& bed,
    GeneticModeSet modes,
    GenotypeMethod method,
    Index chunk_size,
    std::function<void(std::size_t)> observer)
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
    const LocusEncoder encoder{bed_};

    std::vector<EncodingSpec> specs;
    specs.reserve(modes_.size());
    for (GeneticMode mode : modes_.each())
    {
        specs.push_back(encoding_spec_from_method(mode, method_));
    }

    std::vector<Eigen::MatrixXd> grms(
        modes_.size(), Eigen::MatrixXd::Zero(n, n));

    for (Index s = start; s < end; s += chunk_size_)
    {
        const Index e = std::min(s + chunk_size_, end);
        const Index cols = e - s;
        std::vector<Eigen::MatrixXd> z(specs.size(), Eigen::MatrixXd(n, cols));

        // Genotype counts are mode-independent, so each variant is tabulated
        // once and shared across modes; expansion then fills each mode's Z.
        for (Index c = 0; c < cols; ++c)
        {
            const LocusStats stats{encoder.count(s + c)};
            for (std::size_t k = 0; k < specs.size(); ++k)
            {
                const LocusEncoding encoding{
                    encoder.encoding(s + c, stats, specs[k])};
                encoder.expand(s + c, encoding, z[k].col(c));
            }
        }

        // grms[k] += z[k] * z[k]^T, accumulating into the lower triangle.
        for (std::size_t k = 0; k < specs.size(); ++k)
        {
            grms[k].selfadjointView<Eigen::Lower>().rankUpdate(z[k]);
        }

        processed_ += cols;
        notify(observer_, static_cast<std::size_t>(processed_));
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

auto GrmBuilder::build(std::span<const MarkerRange> ranges, const Sink& sink)
    -> void
{
    processed_ = 0;

    for (const auto& range : ranges)
    {
        for (const auto& matrix :
             accumulate(range.label, range.start, range.end))
        {
            sink(matrix);
        }
    }
}

}  // namespace gelex
