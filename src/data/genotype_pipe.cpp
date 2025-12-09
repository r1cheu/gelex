#include "gelex/data/genotype_pipe.h"

#include <filesystem>
#include <utility>

#include "gelex/exception.h"

namespace gelex
{

GenotypePipe::GenotypePipe(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    const std::filesystem::path& output_prefix)
    : bed_pipe_(bed_path, std::move(sample_manager))
{
    auto logger = logging::get();

    auto matrix_path = output_prefix;
    matrix_path += ".bmat";
    auto stats_path = output_prefix;
    stats_path += ".snpstats";

    bool exists = std::filesystem::exists(matrix_path)
                  || std::filesystem::exists(stats_path);

    if (exists)
    {
        logger->error(
            "Output files already exist: [{}] or [{}]",
            matrix_path.string(),
            stats_path.string());
        throw FileExistsException(
            std::format("{}: existing files", output_prefix.string()));
    }

    matrix_writer_ = std::make_unique<detail::BinaryMatrixWriter>(matrix_path);
    stats_writer_ = std::make_unique<detail::SnpStatsWriter>(stats_path);

    num_variants_ = bed_pipe_.num_snps();
    sample_size_ = bed_pipe_.num_samples();
}

GenotypeMap GenotypePipe::finalize()
{
    stats_writer_->write(
        sample_size_, monomorphic_indices_, means_, variances_);

    return GenotypeMap(matrix_writer_->path());
}

}  // namespace gelex
