#include "gelex/data/genotype_pipe.h"

#include <expected>
#include <filesystem>
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
    BedPipe&& bed_pipe,
    detail::BinaryMatrixWriter&& matrix_writer,
    detail::SnpStatsWriter&& stats_writer,
    double monomorphic_threshold)
    : bed_pipe_(std::move(bed_pipe)),
      monomorphic_threshold_(monomorphic_threshold),
      num_variants_(bed_pipe_.num_variants()),
      matrix_writer_(std::move(matrix_writer)),
      stats_writer_(std::move(stats_writer))
{
    means_.reserve(num_variants_);
    stddevs_.reserve(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);
}

GenotypePipe::~GenotypePipe() = default;

auto GenotypePipe::create(
    const std::filesystem::path& bed_prefix,
    bool iid_only,
    double monomorphic_threshold) -> std::expected<GenotypePipe, Error>
{
    auto bed_pipe = BedPipe::create(bed_prefix, iid_only);
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    if (monomorphic_threshold <= 0.0)
    {
        return std::unexpected(
            Error{
                ErrorCode::InvalidData,
                "Monomorphic threshold must be positive"});
    }

    std::filesystem::path matrix_path = bed_prefix;
    matrix_path.replace_extension(".bmat");
    std::filesystem::path stats_path = bed_prefix;
    stats_path.replace_extension(".snpstats");

    auto matrix_writer = detail::BinaryMatrixWriter::create(matrix_path);
    if (!matrix_writer)
    {
        return std::unexpected(matrix_writer.error());
    }
    auto stats_writer = detail::SnpStatsWriter::create(stats_path);
    if (!stats_writer)
    {
        return std::unexpected(stats_writer.error());
    }

    return GenotypePipe{
        std::move(*bed_pipe),
        std::move(*matrix_writer),
        std::move(*stats_writer),
        monomorphic_threshold};
}

auto GenotypePipe::process(
    size_t chunk_size,
    const std::optional<std::unordered_map<std::string, Index>>& process_id_map)
    -> std::expected<void, Error>
{
    auto logger = gelex::logging::get();
    if (process_id_map)
    {
        sample_size_ = static_cast<Index>(process_id_map->size());
    }
    else
    {
        sample_size_ = bed_pipe_.sample_size();
    }

    // Process in chunks
    for (Index start_variant = 0; start_variant < num_variants_;)
    {
        Index end_variant = std::min(
            static_cast<Index>(start_variant + chunk_size), num_variants_);
        auto chunk
            = bed_pipe_.load_chunk(start_variant, end_variant, process_id_map);
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
    if (auto result = matrix_writer_.write(chunk); !result)
    {
        return std::unexpected(result.error());
    }

    return {};
}

auto GenotypePipe::finalize() -> std::expected<void, Error>
{
    // Write final statistics
    if (auto result = stats_writer_.write(
            sample_size_,
            num_variants_,
            monomorphic_indices_.size(),
            monomorphic_indices_,
            means_,
            stddevs_);
        !result)
    {
        return std::unexpected(result.error());
    }
    return {};
}

}  // namespace gelex
