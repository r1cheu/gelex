#include "gelex/data/genotype_pipe.h"

#include <expected>
#include <filesystem>
#include <memory>
#include <optional>
#include <utility>

#include <Eigen/Core>

#include "../src/data/binary_matrix_writer.h"
#include "../src/data/snp_stats_writer.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

using Eigen::Index;

GenotypePipe::GenotypePipe(
    std::unique_ptr<BedPipe> bed_pipe,
    std::filesystem::path matrix_path,
    std::filesystem::path stats_path,
    std::optional<std::unordered_map<std::string, Index>> id_map,
    double monomorphic_threshold)
    : bed_pipe_(std::move(bed_pipe)),
      matrix_path_(std::move(matrix_path)),
      stats_path_(std::move(stats_path)),
      monomorphic_threshold_(monomorphic_threshold),
      id_map_(std::move(id_map))
{
    // Determine sample size based on ID map or original data
    if (id_map_.has_value())
    {
        num_samples_ = static_cast<Index>(id_map_->size());
    }
    else
    {
        num_samples_ = bed_pipe_->sample_size();
    }

    num_variants_ = bed_pipe_->num_variants();

    means_.reserve(num_variants_);
    stddevs_.reserve(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);
}

GenotypePipe::~GenotypePipe() = default;

auto GenotypePipe::create(
    std::unique_ptr<BedPipe> bed_pipe,
    std::filesystem::path matrix_path,
    std::filesystem::path stats_path,
    std::optional<std::unordered_map<std::string, Index>> id_map,
    double monomorphic_threshold) -> std::expected<GenotypePipe, Error>
{
    if (!bed_pipe)
    {
        return std::unexpected(
            Error{ErrorCode::InvalidData, "BedPipe cannot be null"});
    }

    if (monomorphic_threshold <= 0.0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                "Monomorphic threshold must be positive"});
    }

    return GenotypePipe{
        std::move(bed_pipe),
        std::move(matrix_path),
        std::move(stats_path),
        std::move(id_map),
        monomorphic_threshold};
}

auto GenotypePipe::process(
    size_t chunk_size,
    std::optional<std::unordered_map<std::string, Index>> process_id_map)
    -> std::expected<void, Error>
{
    auto logger = gelex::logging::get();

    // Use process-level ID map if provided, otherwise use constructor-level
    auto effective_id_map = process_id_map.has_value()
                                ? std::move(process_id_map)
                                : std::move(id_map_);

    // Update sample size if using process-level ID map
    if (effective_id_map.has_value())
    {
        num_samples_ = static_cast<Index>(effective_id_map->size());
    }

    logger->info(
        "Starting genotype processing: {} samples, {} variants, chunk size {}",
        num_samples_,
        num_variants_,
        chunk_size);

    // Create writers
    matrix_writer_ = std::make_unique<detail::BinaryMatrixWriter>(matrix_path_);
    stats_writer_ = std::make_unique<detail::SnpStatsWriter>(stats_path_);

    if (auto result = matrix_writer_->open(); !result)
    {
        return std::unexpected(result.error());
    }

    // Process in chunks
    for (Index start_variant = 0; start_variant < num_variants_;)
    {
        Index end_variant = std::min(
            static_cast<Index>(start_variant + chunk_size), num_variants_);

        auto chunk = bed_pipe_->load_chunk(
            start_variant, end_variant, effective_id_map);
        if (!chunk)
        {
            return std::unexpected(chunk.error());
        }

        if (auto result = process_chunk(std::move(*chunk), start_variant);
            !result)
        {
            return std::unexpected(result.error());
        }

        start_variant = end_variant;

        if (start_variant % (10 * chunk_size) == 0)
        {
            logger->info("Processed {} variants", start_variant);
        }
    }

    logger->info("Completed processing all {} variants", num_variants_);

    return finalize();
}

auto GenotypePipe::process_chunk(Eigen::MatrixXd&& chunk, size_t global_start)
    -> std::expected<void, Error>
{
    const Index num_variants_in_chunk = chunk.cols();

    for (Index variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        const size_t global_idx = global_start + variant_idx;
        auto variant = chunk.col(variant_idx);

        // Compute statistics in single pass
        const auto n = static_cast<double>(variant.size());
        const double mean = variant.mean();
        const double variance
            = (variant.array() - mean).square().sum() / (n - 1);
        const double std_dev = std::sqrt(variance);

        means_.push_back(mean);
        stddevs_.push_back(std_dev);

        // Check for monomorphic SNP
        if (std_dev < monomorphic_threshold_)
        {
            monomorphic_indices_.push_back(global_idx);
            variant.setZero();  // Standardize to all zeros
        }
        else
        {
            // Standardize: (x - mean) / std_dev
            variant = (variant.array() - mean) / std_dev;
        }
    }

    // Write processed chunk to binary matrix
    if (auto result = matrix_writer_->append_matrix(chunk); !result)
    {
        return std::unexpected(result.error());
    }

    return {};
}

auto GenotypePipe::finalize() -> std::expected<void, Error>
{
    // Write final statistics
    if (auto result = stats_writer_->write_all(
            num_samples_,
            num_variants_,
            monomorphic_indices_.size(),
            monomorphic_indices_,
            means_,
            stddevs_);
        !result)
    {
        return std::unexpected(result.error());
    }

    // Close writers
    if (auto result = matrix_writer_->close(); !result)
    {
        return std::unexpected(result.error());
    }

    if (auto result = stats_writer_->close(); !result)
    {
        return std::unexpected(result.error());
    }

    return {};
}

}  // namespace gelex
