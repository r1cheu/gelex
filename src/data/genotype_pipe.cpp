#include "gelex/data/genotype_pipe.h"

#include <filesystem>
#include <utility>

#include "gelex/exception.h"

namespace gelex
{

GenotypePipe::GenotypePipe(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager,
    const std::filesystem::path& output_prefix,
    bool force_overwrite)
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
        if (!force_overwrite)
        {
            logger->error(
                "Output files already exist: [{}] or [{}]. Use "
                "force_overwrite=true to bypass.",
                matrix_path.string(),
                stats_path.string());
            throw OutputFileExistsException(matrix_path);
        }

        logger->warn(
            "Overwriting existing output files: [{}]", output_prefix.string());
    }

    if (auto parent = matrix_path.parent_path(); !parent.empty())
    {
        std::filesystem::create_directories(parent);
    }

    matrix_writer_ = std::make_unique<detail::BinaryMatrixWriter>(matrix_path);
    stats_writer_ = std::make_unique<detail::SnpStatsWriter>(stats_path);

    num_variants_ = bed_pipe_.num_snps();
    sample_size_ = bed_pipe_.num_samples();
}

GenotypePipe::~GenotypePipe()
{
    try
    {
        wait_for_write();
    }
    catch (...)
    {
    }
}

GenotypeMap GenotypePipe::finalize()
{
    wait_for_write();

    stats_writer_->write(
        sample_size_, monomorphic_indices_, means_, variances_);

    return GenotypeMap(matrix_writer_->path());
}

}  // namespace gelex
