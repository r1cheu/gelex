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

#ifndef GELEX_DATA_GENOTYPE_PIPE_H_
#define GELEX_DATA_GENOTYPE_PIPE_H_

#include <cstdint>
#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Core>
#include "barkeep.h"

#include "../src/data/binary_matrix_writer.h"
#include "../src/data/snp_stats_writer.h"
#include "../src/estimator/bayes/indicator.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_mmap.h"
#include "gelex/data/sample_manager.h"
#include "gelex/data/variant_processor.h"

namespace gelex
{

class GenotypePipe
{
   public:
    GenotypePipe(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager,
        const std::filesystem::path& output_prefix);

    GenotypePipe(const GenotypePipe&) = delete;
    GenotypePipe(GenotypePipe&&) noexcept = default;
    GenotypePipe& operator=(const GenotypePipe&) = delete;
    GenotypePipe& operator=(GenotypePipe&&) noexcept = default;
    ~GenotypePipe() = default;

    template <VariantProcessor Processor = StandardizingProcessor>
    auto process(size_t chunk_size = 10000) -> GenotypeMap
    {
        int64_t current_processed_snps = 0;
        Processor processor;
        auto pbar = detail::create_genotype_process_bar<Processor>(
            current_processed_snps, num_variants_);
        means_.resize(num_variants_);
        variances_.resize(num_variants_);
        monomorphic_indices_.clear();
        monomorphic_indices_.reserve(num_variants_ / 100);

        pbar->show();
        for (int64_t start_variant = 0; start_variant < num_variants_;)
        {
            int64_t end_variant = std::min(
                static_cast<int64_t>(start_variant + chunk_size),
                num_variants_);

            auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);

            process_chunk<Processor>(chunk, start_variant, processor);
            matrix_writer_->write(chunk);
            current_processed_snps += (end_variant - start_variant);
            start_variant = end_variant;
        }
        pbar->done();

        return finalize();
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
    template <typename Processor>
    void process_chunk(
        Eigen::MatrixXd& chunk,
        size_t global_start,
        Processor& processor)
    {
        const int64_t num_variants_in_chunk = chunk.cols();

        for (int64_t i = 0; i < num_variants_in_chunk; ++i)
        {
            auto variant = chunk.col(i);

            VariantStats stats = processor.process_variant(variant);

            size_t global_idx = global_start + i;
            means_[global_idx] = stats.mean;
            variances_[global_idx] = stats.stddev;

            if (stats.is_monomorphic)
            {
                monomorphic_indices_.push_back(
                    static_cast<int64_t>(global_idx));
            }
        }
    }

    GenotypeMap finalize();

    BedPipe bed_pipe_;
    int64_t sample_size_{};
    int64_t num_variants_{};

    std::vector<double> means_;
    std::vector<double> variances_;
    std::vector<int64_t> monomorphic_indices_;

    std::unique_ptr<detail::BinaryMatrixWriter> matrix_writer_;
    std::unique_ptr<detail::SnpStatsWriter> stats_writer_;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_PIPE_H_
