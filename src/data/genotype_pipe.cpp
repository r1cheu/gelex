#include "gelex/data/genotype_pipe.h"

#include <cstdint>
#include <expected>
#include <filesystem>
#include <utility>

#include <Eigen/Core>

#include "binary_matrix_writer.h"
#include "estimator/bayes/indicator.h"
#include "gelex/barkeep.h"
#include "gelex/error.h"
#include "gelex/logger.h"
#include "snp_stats_writer.h"

namespace gelex
{

namespace bk = barkeep;
using Eigen::Index;

GenotypePipe::GenotypePipe(
    BedPipe&& bed_pipe,
    detail::BinaryMatrixWriter&& matrix_writer,
    detail::SnpStatsWriter&& stats_writer,
    bool dominant)
    : bed_pipe_(std::move(bed_pipe)),
      monomorphic_threshold_(1e-6),
      dominant_(dominant),
      matrix_writer_(std::move(matrix_writer)),
      stats_writer_(std::move(stats_writer))
{
    num_variants_ = bed_pipe_.num_variants();  // NOLINT
    sample_size_ = bed_pipe_.sample_size();    // NOLINT

    means_.reserve(num_variants_);
    stddevs_.reserve(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);
}

GenotypePipe::~GenotypePipe() = default;

auto GenotypePipe::create(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    const std::filesystem::path& output_prefix,
    bool dominant) -> std::expected<GenotypePipe, Error>
{
    auto bed_pipe = BedPipe::create(bed_path, std::move(sample_manager));
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    std::filesystem::path matrix_path = output_prefix;
    matrix_path.replace_extension(dominant ? ".dom.bmat" : ".add.bmat");
    std::filesystem::path stats_path = output_prefix;
    stats_path.replace_extension(dominant ? ".dom.snpstats" : ".add.snpstats");

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
        dominant};
}

auto GenotypePipe::process(size_t chunk_size) -> std::expected<void, Error>
{
    auto logger = gelex::logging::get();
    global_snp_idx_ = 0;

    auto pbar = bk::ProgressBar(
        &global_snp_idx_,
        {.total = static_cast<uint64_t>(num_variants_),
         .format = "{bar} {value}/{total} ",
         .style = detail::BAR_STYLE});

    for (int64_t start_variant = 0; start_variant < num_variants_;)
    {
        int64_t end_variant = std::min(
            static_cast<int64_t>(start_variant + chunk_size), num_variants_);

        auto chunk
            = bed_pipe_.load_chunk(start_variant, end_variant, dominant_);
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
    pbar->done();

    return finalize();
}

auto GenotypePipe::process_chunk(Eigen::MatrixXd&& chunk, size_t global_start)
    -> std::expected<void, Error>
{
    Eigen::MatrixXd matrix = std::move(chunk);
    const int64_t num_variants_in_chunk = matrix.cols();

    for (int64_t variant_idx = 0; variant_idx < num_variants_in_chunk;
         ++variant_idx)
    {
        global_snp_idx_++;
        const size_t global_idx = global_start + variant_idx;
        auto variant = matrix.col(variant_idx);

        const auto n = static_cast<double>(variant.size());
        const double mean = variant.mean();
        const double variance
            = (variant.array() - mean).square().sum() / (n - 1);
        const double std_dev = std::sqrt(variance);

        means_.push_back(mean);
        stddevs_.push_back(std_dev);

        if (std_dev < monomorphic_threshold_)
        {
            monomorphic_indices_.push_back(static_cast<int64_t>(global_idx));
        }
        else
        {
            variant = (variant.array() - mean) / std_dev;
        }
    }

    if (auto result = matrix_writer_.write(matrix); !result)
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
            static_cast<int64_t>(monomorphic_indices_.size()),
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
