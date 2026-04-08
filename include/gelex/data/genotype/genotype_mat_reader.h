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

#ifndef GELEX_DATA_GENOTYPE_MAT_READER_H_
#define GELEX_DATA_GENOTYPE_MAT_READER_H_

#include <algorithm>
#include <filesystem>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/data/genotype/genotype_matrix.h"
#include "gelex/data/genotype/genotype_processor.h"
#include "gelex/infra/logging/data_pipe_event.h"
#include "gelex/infra/logging/notify.h"
#include "gelex/types/genotype_process_method.h"

namespace gelex
{

class GenotypeMatReader
{
   public:
    explicit GenotypeMatReader(
        const std::filesystem::path& bed_path,
        const df::Index<std::string>& sample_index,
        DataPipeObserver observer = {});

    GenotypeMatReader(const GenotypeMatReader&) = delete;
    GenotypeMatReader& operator=(const GenotypeMatReader&) = delete;
    GenotypeMatReader(GenotypeMatReader&&) noexcept = default;
    GenotypeMatReader& operator=(GenotypeMatReader&&) noexcept = default;
    ~GenotypeMatReader() = default;

    template <GeneticKind GT>
    auto process(GenotypeProcessMethod method, size_t chunk_size = 10000)
        -> GenotypeMatrix
    {
        global_snp_idx_ = 0;
        auto fn = get_genotype_process_method<GT>(method);
        means_.resize(num_variants_);
        stddevs_.resize(num_variants_);
        monomorphic_indices_.reserve(num_variants_ / 100);

        for (int64_t start_variant = 0; start_variant < num_variants_;)
        {
            int64_t end_variant = std::min(
                static_cast<int64_t>(start_variant + chunk_size),
                num_variants_);
            auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);
            process_chunk(chunk, start_variant, fn);
            global_snp_idx_ += chunk.cols();
            notify(
                observer_,
                GenotypeProgressEvent{
                    static_cast<size_t>(global_snp_idx_),
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
    void process_chunk(
        Eigen::MatrixXd& chunk,
        Eigen::Index global_start,
        LocusStatistic (*fn)(Eigen::Ref<Eigen::VectorXd>));

    GenotypeMatrix finalize();

    BedPipe bed_pipe_;
    DataPipeObserver observer_;

    int64_t sample_size_{};
    int64_t num_variants_{};

    int64_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> stddevs_;
    std::vector<int64_t> monomorphic_indices_;

    Eigen::MatrixXd data_matrix_;
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_MAT_READER_H_
