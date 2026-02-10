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

#ifndef GELEX_DATA_GENOTYPE_LOADER_H_
#define GELEX_DATA_GENOTYPE_LOADER_H_

#include <algorithm>
#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_processor.h"
#include "gelex/data/sample_manager.h"
#include "gelex/detail/indicator.h"

namespace gelex
{

class GenotypeLoader
{
   public:
    explicit GenotypeLoader(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager);

    GenotypeLoader(const GenotypeLoader&) = delete;
    GenotypeLoader& operator=(const GenotypeLoader&) = delete;
    GenotypeLoader(GenotypeLoader&&) noexcept = default;
    GenotypeLoader& operator=(GenotypeLoader&&) noexcept = default;
    ~GenotypeLoader() = default;

    template <
        GenotypeProcessor Processor = AdditiveProcessor<StandardizeMethod>>
    GenotypeMatrix process(size_t chunk_size = 10000);

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
        Eigen::Index global_start,
        Processor& processor);

    GenotypeMatrix finalize();

    BedPipe bed_pipe_;

    int64_t sample_size_{};
    int64_t num_variants_{};

    int64_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> stddevs_;
    std::vector<int64_t> monomorphic_indices_;

    Eigen::MatrixXd data_matrix_;
};

template <GenotypeProcessor Processor>
GenotypeMatrix GenotypeLoader::process(size_t chunk_size)
{
    global_snp_idx_ = 0;

    Processor processor;
    auto pbar = detail::create_genotype_process_bar<Processor>(
        global_snp_idx_, num_variants_);
    pbar->show();
    means_.resize(num_variants_);
    stddevs_.resize(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);

    for (int64_t start_variant = 0; start_variant < num_variants_;)
    {
        int64_t end_variant = std::min(
            static_cast<int64_t>(start_variant + chunk_size), num_variants_);

        auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);
        process_chunk(chunk, start_variant, processor);

        global_snp_idx_ += chunk.cols();
        start_variant = end_variant;
    }
    pbar->done();
    return finalize();
}

template <typename Processor>
void GenotypeLoader::process_chunk(
    Eigen::MatrixXd& chunk,
    Eigen::Index global_start,
    Processor& processor)
{
    const Eigen::Index num_variants_in_chunk = chunk.cols();

#pragma omp parallel for schedule(static)
    for (Eigen::Index i = 0; i < num_variants_in_chunk; ++i)
    {
        auto variant = chunk.col(i);
        VariantStats stats = processor.process_variant(variant);

        Eigen::Index global_idx = global_start + i;
        means_[global_idx] = stats.mean;
        stddevs_[global_idx] = stats.stddev;

        if (stats.is_monomorphic)
        {
#pragma omp critical(genotype_loader_mono)
            {
                monomorphic_indices_.push_back(
                    static_cast<int64_t>(global_idx));
            }
        }
    }

    data_matrix_.middleCols(global_start, num_variants_in_chunk) = chunk;
}

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_LOADER_H_
