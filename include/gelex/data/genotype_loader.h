#pragma once

#include <expected>
#include <filesystem>
#include <memory>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/genotype_matrix.h"
#include "gelex/data/genotype_pipe.h"  // For VariantStats and VariantProcessor
#include "gelex/data/sample_manager.h"
#include "gelex/error.h"

namespace gelex
{

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
     * @brief Load genotype data from BED file into memory
     *
     * @param bed_path Path to PLINK BED file
     * @param sample_manager Sample manager for filtering/ordering samples
     * @param dominant Load dominant encoding (default: false for additive)
     * @param chunk_size Number of variants to process at once
     * @return GenotypeMatrix with processed data
     */
    static auto load_from_bed(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager,
        bool dominant = false,
        size_t chunk_size = 10000) -> std::expected<GenotypeMatrix, Error>;

    /**
     * @brief Load with custom variant processor
     *
     * @tparam Processor Variant processing strategy
     * @param bed_path Path to PLINK BED file
     * @param sample_manager Sample manager for filtering/ordering samples
     * @param dominant Load dominant encoding
     * @param chunk_size Number of variants to process at once
     * @return GenotypeMatrix with processed data
     */
    template <VariantProcessor Processor>
    static auto load_with_processor(
        const std::filesystem::path& bed_path,
        std::shared_ptr<SampleManager> sample_manager,
        bool dominant = false,
        size_t chunk_size = 10000) -> std::expected<GenotypeMatrix, Error>;

   private:
    // Internal implementation for template method
    template <VariantProcessor Processor>
    static auto load_impl(
        BedPipe& bed_pipe,
        bool dominant,
        size_t chunk_size,
        Processor& processor) -> std::expected<GenotypeMatrix, Error>;
};

// Template implementation
template <VariantProcessor Processor>
auto GenotypeLoader::load_with_processor(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    bool dominant,
    size_t chunk_size) -> std::expected<GenotypeMatrix, Error>
{
    auto bed_pipe = BedPipe::create(bed_path, std::move(sample_manager));
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    Processor processor;
    return load_impl(*bed_pipe, dominant, chunk_size, processor);
}

template <VariantProcessor Processor>
auto GenotypeLoader::load_impl(
    BedPipe& bed_pipe,
    bool dominant,
    size_t chunk_size,
    Processor& processor) -> std::expected<GenotypeMatrix, Error>
{
    const int64_t num_samples = bed_pipe.sample_size();
    const int64_t num_variants = bed_pipe.num_variants();

    // Pre-allocate matrix for all genotype data
    Eigen::MatrixXd data(num_samples, num_variants);

    // Storage for statistics
    std::vector<double> means;
    std::vector<double> variances;
    std::vector<int64_t> monomorphic_indices;

    means.reserve(static_cast<size_t>(num_variants));
    variances.reserve(static_cast<size_t>(num_variants));
    monomorphic_indices.reserve(static_cast<size_t>(num_variants) / 100);

    // Process variants in chunks
    int64_t global_variant_idx = 0;

    for (int64_t start_variant = 0; start_variant < num_variants;)
    {
        const int64_t end_variant = std::min(
            start_variant + static_cast<int64_t>(chunk_size), num_variants);

        auto chunk = bed_pipe.load_chunk(start_variant, end_variant, dominant);
        if (!chunk)
        {
            return std::unexpected(chunk.error());
        }

        const int64_t num_variants_in_chunk = chunk->cols();

        // Process each variant in the chunk
        for (int64_t variant_idx = 0; variant_idx < num_variants_in_chunk;
             ++variant_idx)
        {
            auto variant = chunk->col(variant_idx);

            VariantStats stats = processor.process_variant(variant);

            data.col(global_variant_idx) = variant;

            // Store statistics
            means.push_back(stats.mean);
            variances.push_back(stats.variance);

            if (stats.is_monomorphic)
            {
                monomorphic_indices.push_back(global_variant_idx);
            }

            ++global_variant_idx;
        }

        start_variant = end_variant;
    }

    // Convert statistics to Eigen vectors
    Eigen::VectorXd mean_vec
        = Eigen::Map<Eigen::VectorXd>(means.data(), means.size());
    Eigen::VectorXd variance_vec
        = Eigen::Map<Eigen::VectorXd>(variances.data(), variances.size());

    // Convert monomorphic indices to set
    std::unordered_set<int64_t> mono_set(
        monomorphic_indices.begin(), monomorphic_indices.end());

    return GenotypeMatrix::create(
        std::move(data),
        std::move(mono_set),
        std::move(mean_vec),
        std::move(variance_vec));
}

}  // namespace gelex
