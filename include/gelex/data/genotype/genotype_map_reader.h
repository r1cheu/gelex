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

#ifndef GELEX_DATA_GENOTYPE_MAP_READER_H_
#define GELEX_DATA_GENOTYPE_MAP_READER_H_

#include <fmt/format.h>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/genotype/genotype_mmap.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/io/binary_writer.h"
#include "gelex/types/genotype_process_method.h"

namespace gelex
{

class GenotypeMapReader
{
   public:
    GenotypeMapReader(
        const std::filesystem::path& bed_path,
        const df::Index<std::string>& sample_index,
        const std::filesystem::path& output_prefix,
        DataPipeObserver observer = {});

    template <GeneticMode GT>
    auto process(GenotypeProcessMethod method, size_t chunk_size = 10000)
        -> GenotypeMap
    {
        constexpr auto effect = EffectType::from_genetic(GT);
        const auto geno_path = fmt::format("{}/genotype", effect);
        genotype_handle_
            = writer_->reserve<double>(geno_path, sample_size_, num_variants_);

        int64_t current_processed_snps = 0;
        auto fn = get_genotype_process_method<GT>(method);
        means_.resize(num_variants_);
        variances_.resize(num_variants_);
        monomorphic_indices_.clear();
        monomorphic_indices_.reserve(num_variants_ / 100);

        for (int64_t start_variant = 0; start_variant < num_variants_;)
        {
            int64_t end_variant = std::min(
                static_cast<int64_t>(start_variant + chunk_size),
                num_variants_);

            auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);
            process_chunk(chunk, start_variant, fn);
            writer_->write(genotype_handle_, chunk);
            current_processed_snps += (end_variant - start_variant);

            notify(
                observer_,
                GenotypeProgressEvent{
                    static_cast<size_t>(current_processed_snps),
                    static_cast<size_t>(num_variants_),
                    false});

            start_variant = end_variant;
        }
        notify(
            observer_,
            GenotypeProgressEvent{
                static_cast<size_t>(num_variants_),
                static_cast<size_t>(num_variants_),
                true});

        return finalize<GT>();
    }

    [[nodiscard]] Eigen::Index num_samples() const noexcept
    {
        return sample_size_;
    }
    [[nodiscard]] Eigen::Index num_variants() const noexcept
    {
        return num_variants_;
    }

   private:
    void process_chunk(
        Eigen::MatrixXd& chunk,
        size_t global_start,
        LocusStatistic (*fn)(Eigen::Ref<Eigen::VectorXd>));

    template <GeneticMode GT>
    auto finalize() -> GenotypeMap
    {
        constexpr auto effect = EffectType::from_genetic(GT);

        auto stats_handle = writer_->reserve<double>(
            fmt::format("{}/loci_stats", effect), num_variants_, 2);
        writer_->write(stats_handle, means_);
        writer_->write(stats_handle, variances_);

        if (!monomorphic_indices_.empty())
        {
            auto mono_handle = writer_->reserve<int64_t>(
                fmt::format("{}/mono_indices", effect),
                monomorphic_indices_.size(),
                1);
            writer_->write(mono_handle, monomorphic_indices_);
        }

        writer_->finalize();

        auto gbin_path = output_prefix_;
        gbin_path += ".gbin";
        return GenotypeMap(gbin_path, GT);
    }

    BedPipe bed_pipe_;
    DataPipeObserver observer_;
    int64_t sample_size_{};
    int64_t num_variants_{};

    std::vector<double> means_;
    std::vector<double> variances_;
    std::vector<int64_t> monomorphic_indices_;

    std::unique_ptr<detail::BinaryWriter> writer_;
    detail::SectionHandle<double> genotype_handle_{};
    std::filesystem::path output_prefix_;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_MAP_READER_H_
