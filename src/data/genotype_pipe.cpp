#include "gelex/data/genotype_pipe.h"

#include <cassert>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <utility>

#include <Eigen/Core>

#include "binary_matrix_writer.h"
#include "gelex/data/genotype_mmap.h"
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

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / std_dev;
    }

    return stats;
}

VariantStats RawProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;
    const auto n = static_cast<double>(variant.size());

    stats.mean = variant.mean();
    stats.variance = (variant.array() - stats.mean).square().sum() / (n - 1);

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    return stats;
}

VariantStats HardWenbergProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;

    stats.mean = variant.mean();
    stats.variance = stats.mean * (1 - 0.5 * stats.mean);
    const double std_dev = std::sqrt(stats.variance);

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / std_dev;
    }

    return stats;
}

VariantStats DominantStandardizingProcessor::process_variant(
    Ref<VectorXd> variant)
{
    VariantStats stats;
    const auto n = static_cast<double>(variant.size());

    variant = (variant.array() == 2).select(0, variant.array());

    stats.mean = variant.mean();
    stats.variance = (variant.array() - stats.mean).square().sum() / (n - 1);
    const double std_dev = std::sqrt(stats.variance);

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / std_dev;
    }

    return stats;
}

VariantStats DominantRawProcessor::process_variant(Ref<VectorXd> variant)
{
    VariantStats stats;
    const auto n = static_cast<double>(variant.size());

    variant = (variant.array() == 2).select(0, variant.array());

    stats.mean = variant.mean();
    stats.variance = (variant.array() - stats.mean).square().sum() / (n - 1);

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    return stats;
}

VariantStats DominantOrthogonalHWEProcessor::process_variant(
    Ref<VectorXd> variant)
{
    VariantStats stats;
    const double p_freq = variant.mean() / 2;

    stats.mean = 2 * p_freq * p_freq;
    const double stddev = 2 * p_freq * (1 - p_freq);
    stats.variance = stddev * stddev;

    const double one_alt_encode = 2 * p_freq;
    const double two_alt_encode = (4 * p_freq) - 2;
    auto op = [&one_alt_encode, &two_alt_encode](double x) -> double
    {
        if (x == 1.0)
        {
            return one_alt_encode;
        }
        if (x == 2.0)
        {
            return two_alt_encode;
        }
        return x;
    };
    variant = variant.unaryExpr(op);

    stats.is_monomorphic
        = (stats.variance - 0) < std::numeric_limits<double>::epsilon();

    if (!stats.is_monomorphic)
    {
        variant = (variant.array() - stats.mean) / stddev;
    }

    return stats;
}

GenotypePipe::GenotypePipe(
    BedPipe&& bed_pipe,
    detail::BinaryMatrixWriter&& matrix_writer,
    detail::SnpStatsWriter&& stats_writer)
    : bed_pipe_(std::move(bed_pipe)),
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
    const std::filesystem::path& output_prefix)
    -> std::expected<GenotypePipe, Error>
{
    auto logger = logging::get();
    auto bed_pipe = BedPipe::create(bed_path, std::move(sample_manager));
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    std::filesystem::path matrix_path = output_prefix;
    matrix_path += ".bmat";
    std::filesystem::path stats_path = output_prefix;
    stats_path += ".snpstats";

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
        std::move(*stats_writer)};
}

auto GenotypePipe::finalize() -> std::expected<GenotypeMap, Error>
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

    return GenotypeMap::create(matrix_writer_.path());
}

}  // namespace gelex
