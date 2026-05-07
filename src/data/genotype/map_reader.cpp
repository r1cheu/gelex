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

#include "gelex/data/genotype/map_reader.h"

#include <cstdint>
#include <filesystem>
#include <utility>

#include <Eigen/Core>

#include "gelex/data/genotype/process_method.h"
#include "gelex/exception.h"
#include "gelex/infra/logger.h"

namespace gelex::genotype
{

void GenotypeMapReader::process_chunk(
    Eigen::MatrixXd& chunk,
    size_t global_start,
    gelex::LocusStatistic (*fn)(Eigen::Ref<Eigen::VectorXd>))
{
    const int64_t num_variants_in_chunk = chunk.cols();

    for (int64_t i = 0; i < num_variants_in_chunk; ++i)
    {
        auto variant = chunk.col(i);
        gelex::LocusStatistic stats = fn(variant);

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
    const dataframe::Index<std::string>& sample_index,
    const std::filesystem::path& output_prefix,
    gelex::GenoObserver observer)
    : bed_pipe_(bed_path, sample_index),
      observer_(std::move(observer)),
      output_prefix_(output_prefix)
{
    auto gbin_path = output_prefix;
    gbin_path += ".gbin";

    if (std::filesystem::exists(gbin_path))
    {
        auto logger = gelex::logging::get();
        logger->error("Output file already exists: [{}]", gbin_path.string());
        throw gelex::GelexException(
            fmt::format("{}: existing file", gbin_path.string()));
    }

    num_variants_ = bed_pipe_.num_snps();
    sample_size_ = bed_pipe_.num_samples();
}

}  // namespace gelex::genotype
