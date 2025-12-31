#pragma once

#include <algorithm>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"
#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/sample_manager.h"
#include "gelex/data/variant_processor.h"
#include "gelex/logger.h"

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

    template <VariantProcessor Processor = StandardizingProcessor>
    GenotypeMatrix process(
        size_t chunk_size = 10000,
        std::function<void(size_t processed, size_t total)> progress_callback
        = nullptr);

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

    // 成员变量
    BedPipe bed_pipe_;  // 必须首先初始化

    int64_t sample_size_{};
    int64_t num_variants_{};

    size_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> stddevs_;
    std::vector<int64_t> monomorphic_indices_;

    Eigen::MatrixXd data_matrix_;
};

template <VariantProcessor Processor>
GenotypeMatrix GenotypeLoader::process(
    size_t chunk_size,
    std::function<void(size_t, size_t)> progress_callback)
{
    global_snp_idx_ = 0;

    Processor processor;

    means_.resize(num_variants_);
    stddevs_.resize(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);

    // Initial progress update
    if (progress_callback)
    {
        progress_callback(0, num_variants_);
    }

    for (int64_t start_variant = 0; start_variant < num_variants_;)
    {
        int64_t end_variant = std::min(
            static_cast<int64_t>(start_variant + chunk_size), num_variants_);

        auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);
        process_chunk(chunk, start_variant, processor);

        global_snp_idx_ += chunk.cols();
        start_variant = end_variant;

        if (progress_callback)
        {
            progress_callback(global_snp_idx_, num_variants_);
        }
    }

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
