#include "gelex/data/genotype_loader.h"

#include <format>
#include <new>  // for std::bad_alloc
#include "Eigen/Core"

namespace gelex
{

GenotypeLoader::GenotypeLoader(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager)
    : bed_pipe_(bed_path, std::move(sample_manager))
{
    num_variants_ = bed_pipe_.num_snps();    // NOLINT
    sample_size_ = bed_pipe_.num_samples();  // NOLINT

    try
    {
        data_matrix_.resize(sample_size_, num_variants_);
    }
    catch (const std::bad_alloc&)
    {
        throw std::runtime_error(
            std::format(
                "Memory allocation failed for Genotype Matrix ({} x {}). "
                "Requires approx {:.2f} GB RAM.",
                sample_size_,
                num_variants_,
                (double)sample_size_ * num_variants_ * sizeof(double) / 1024.0
                    / 1024.0 / 1024.0));
    }
}

GenotypeMatrix GenotypeLoader::finalize()
{
    // 将 std::vector 映射为 Eigen::VectorXd
    Eigen::VectorXd mean_vec = Eigen::Map<Eigen::VectorXd>(
        means_.data(), static_cast<Eigen::Index>(means_.size()));
    Eigen::VectorXd stddev_vec = Eigen::Map<Eigen::VectorXd>(
        stddevs_.data(), static_cast<Eigen::Index>(stddevs_.size()));

    return GenotypeMatrix(
        std::move(data_matrix_),
        std::move(monomorphic_indices_),
        std::move(mean_vec),
        std::move(stddev_vec));
}

}  // namespace gelex
