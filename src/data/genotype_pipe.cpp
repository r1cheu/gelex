#include "gelex/data/genotype_pipe.h"

#include <cstdint>
#include <expected>
#include <filesystem>
#include <utility>

#include <Eigen/Core>

#include "binary_matrix_writer.h"
#include "gelex/error.h"
#include "gelex/logger.h"
#include "snp_stats_writer.h"

namespace gelex
{

namespace bk = barkeep;
using Eigen::Index;
using Eigen::Ref;
using Eigen::VectorXd;

VariantStats StandardizingProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;
    const auto n = static_cast<double>(variant.size());

    stats.mean = variant.mean();
    stats.variance = (variant.array() - stats.mean).square().sum() / (n - 1);
    const double std_dev = std::sqrt(stats.variance);

    stats.is_monomorphic = (stats.variance - 0) < 1e-8;

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / std_dev;
    }

    return stats;
}

VariantStats NonStandardizingProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;
    const auto n = static_cast<double>(variant.size());

    stats.mean = variant.mean();
    stats.variance = (variant.array() - stats.mean).square().sum() / (n - 1);

    stats.is_monomorphic = (stats.variance - 0) < 1e-8;

    return stats;
}

VariantStats HardWenbergProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;

    stats.mean = variant.mean();
    stats.variance = stats.mean * (1 - 0.5 * stats.mean);
    const double std_dev = std::sqrt(stats.variance);

    stats.is_monomorphic = (stats.variance - 0) < 1e-8;

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / std_dev;
    }

    return stats;
}

GenotypePipe::GenotypePipe(
    BedPipe&& bed_pipe,
    detail::BinaryMatrixWriter&& matrix_writer,
    detail::SnpStatsWriter&& stats_writer,
    bool dominant)
    : bed_pipe_(std::move(bed_pipe)),
      dominant_(dominant),
      matrix_writer_(std::move(matrix_writer)),
      stats_writer_(std::move(stats_writer))
{
    num_variants_ = bed_pipe_.num_variants();  // NOLINT
    sample_size_ = bed_pipe_.sample_size();    // NOLINT

    means_.reserve(num_variants_);
    variances_.reserve(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);
}

auto GenotypePipe::create(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    const std::filesystem::path& output_prefix,
    bool dominant) -> std::expected<GenotypePipe, Error>
{
    auto logger = logging::get();
    auto bed_pipe = BedPipe::create(bed_path, std::move(sample_manager));
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    std::filesystem::path matrix_path = output_prefix;
    matrix_path.replace_extension(dominant ? ".dom.bmat" : ".add.bmat");
    std::filesystem::path stats_path = output_prefix;
    stats_path.replace_extension(dominant ? ".dom.snpstats" : ".add.snpstats");

    if (std::filesystem::exists(matrix_path)
        || std::filesystem::exists(stats_path))
    {
        logger->warn(
            "[{}] or [{}] file exists. "
            "Make sure you are using the same dataset.",
            matrix_path.string(),
            stats_path.string());

        return std::unexpected(
            Error{
                .code = ErrorCode::OutputFileExists,
                .message = "Output file already exists"});
    }

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

auto GenotypePipe::finalize() -> std::expected<void, Error>
{
    // Write final statistics
    if (auto result = stats_writer_.write(
            sample_size_,
            num_variants_,
            static_cast<int64_t>(monomorphic_indices_.size()),
            monomorphic_indices_,
            means_,
            variances_);
        !result)
    {
        return std::unexpected(result.error());
    }
    return {};
}

}  // namespace gelex
