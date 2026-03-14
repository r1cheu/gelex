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

#include "gelex/data/genotype/genotype_map_reader.h"

#include <cstdint>
#include <filesystem>
#include <utility>

#include <Eigen/Core>

#include "gelex/exception.h"
#include "gelex/infra/logger.h"
#include "gelex/types/genotype_process_method.h"

namespace gelex
{

void GenotypeMapReader::process_chunk(
    Eigen::MatrixXd& chunk,
    size_t global_start,
    LocusStatistic (*fn)(Eigen::Ref<Eigen::VectorXd>))
{
    const int64_t num_variants_in_chunk = chunk.cols();

    for (int64_t i = 0; i < num_variants_in_chunk; ++i)
    {
        auto variant = chunk.col(i);
        LocusStatistic stats = fn(variant);

        size_t global_idx = global_start + i;
        means_[global_idx] = stats.mean;
        variances_[global_idx] = stats.stddev;

        if (stats.is_monomorphic)
        {
            monomorphic_indices_.push_back(static_cast<int64_t>(global_idx));
        }
    }
}

GenotypeMapReader::GenotypeMapReader(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    const std::filesystem::path& output_prefix,
    DataPipeObserver observer)
    : bed_pipe_(bed_path, std::move(sample_manager)),
      observer_(std::move(observer)),
      output_prefix_(output_prefix)
{
    auto logger = logging::get();

    auto gbin_path = output_prefix;
    gbin_path += ".gbin";

    if (std::filesystem::exists(gbin_path))
    {
        logger->error("Output file already exists: [{}]", gbin_path.string());
        throw FileExistsException(
            std::format("{}: existing file", gbin_path.string()));
    }

    num_variants_ = bed_pipe_.num_snps();
    sample_size_ = bed_pipe_.num_samples();

    writer_ = std::make_unique<detail::BinaryWriter>(gbin_path.string());
}

}  // namespace gelex
