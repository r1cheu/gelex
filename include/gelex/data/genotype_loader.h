#pragma once

#include <expected>
#include <filesystem>
#include <memory>
#include <vector>

#include <gelex/barkeep.h>
#include <Eigen/Core>

#include "../src/estimator/bayes/indicator.h"  // For BAR_STYLE
#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_pipe.h"  // For VariantStats and VariantProcessor
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

namespace bk = barkeep;

/**
 * @class GenotypeLoader
 * @brief Load genotype data directly into memory without intermediate files
 *
 * Alternative to GenotypePipe that loads entire genotype dataset into memory
 * as a GenotypeMatrix instead of writing to disk and memory-mapping.
 * Suitable for smaller datasets or when minimizing file I/O is desired.
 */
class GenotypeLoader
{
   public:
    /**
     * @brief Create a GenotypeLoader instance
     *
     * @param bed_path Path to PLINK BED file
     * @param sample_manager Sample manager for filtering/ordering samples
     * @param dominant Load dominant encoding (default: false for additive)
     * @return GenotypeLoader instance
     */
    static auto create(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager)
        -> std::expected<GenotypeLoader, Error>;

    GenotypeLoader(const GenotypeLoader&) = delete;
    GenotypeLoader(GenotypeLoader&&) noexcept = default;
    GenotypeLoader& operator=(const GenotypeLoader&) = delete;
    GenotypeLoader& operator=(GenotypeLoader&&) noexcept = default;
    ~GenotypeLoader() = default;

    /**
     * @brief Process genotype data with progress tracking
     *
     * @tparam Processor Variant processing strategy
     * @param chunk_size Number of variants to process at once
     * @return Success or error
     */
    template <VariantProcessor Processor = StandardizingProcessor>
    auto process(size_t chunk_size = 10000)
        -> std::expected<GenotypeMatrix, Error>;

    // Dimension accessors
    Eigen::Index num_samples() const noexcept { return sample_size_; }
    Eigen::Index num_variants() const noexcept { return num_variants_; }
    size_t num_processed_variants() const noexcept { return means_.size(); }

   private:
    explicit GenotypeLoader(BedPipe&& bed_pipe);

    template <VariantProcessor Processor>
    auto process_chunk(
        Eigen::MatrixXd&& chunk,
        Eigen::Index global_start,
        Processor& processor) -> std::expected<void, Error>;

    auto finalize() -> std::expected<GenotypeMatrix, Error>;

    BedPipe bed_pipe_;

    int64_t sample_size_{};
    int64_t num_variants_{};

    size_t global_snp_idx_{};

    std::vector<double> means_;
    std::vector<double> variances_;
    std::vector<int64_t> monomorphic_indices_;

    // Internal storage for building the final matrix
    Eigen::MatrixXd data_matrix_;
};

// Template implementation for process method
template <VariantProcessor Processor>
auto GenotypeLoader::process(size_t chunk_size)
    -> std::expected<GenotypeMatrix, Error>
{
    auto logger = gelex::logging::get();
    global_snp_idx_ = 0;

    logger->info("");
    logger->info("Loading genotype data into memory...");

    auto pbar = bk::ProgressBar(
        &global_snp_idx_,
        {.total = static_cast<uint64_t>(num_variants_),
         .format = "{bar} {value}/{total} ",
         .style = detail::BAR_STYLE});

    Processor processor;

    for (int64_t start_variant = 0; start_variant < num_variants_;)
    {
        int64_t end_variant = std::min(
            static_cast<int64_t>(start_variant + chunk_size), num_variants_);

        auto chunk = bed_pipe_.load_chunk(start_variant, end_variant);
        if (!chunk)
        {
            return std::unexpected(chunk.error());
        }

        if (auto result
            = process_chunk(std::move(*chunk), start_variant, processor);
            !result)
        {
            return std::unexpected(result.error());
        }

        start_variant = end_variant;
    }
    pbar->done();

    return finalize();
}

template <VariantProcessor Processor>
auto GenotypeLoader::process_chunk(
    Eigen::MatrixXd&& chunk,
    Eigen::Index global_start,
    Processor& processor) -> std::expected<void, Error>
{
    Eigen::MatrixXd matrix = std::move(chunk);
    const int64_t num_variants_in_chunk = matrix.cols();

    for (Eigen::Index variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        global_snp_idx_++;
        const Eigen::Index global_idx = global_start + variant_idx;
        auto variant = matrix.col(variant_idx);

        VariantStats stats = processor.process_variant(variant);

        // Store processed variant in data matrix
        data_matrix_.col(global_idx) = variant;

        // Store statistics
        means_.push_back(stats.mean);
        variances_.push_back(stats.variance);

        if (stats.is_monomorphic)
        {
            monomorphic_indices_.push_back(static_cast<int64_t>(global_idx));
        }
    }

    return {};
}
}  // namespace gelex
