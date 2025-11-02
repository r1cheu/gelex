#include "gelex/data/genotype_loader.h"

#include <gelex/barkeep.h>

#include "gelex/logger.h"

namespace gelex
{

namespace bk = barkeep;

GenotypeLoader::GenotypeLoader(BedPipe&& bed_pipe)
    : bed_pipe_(std::move(bed_pipe))
{
    num_variants_ = bed_pipe_.num_variants();
    sample_size_ = bed_pipe_.sample_size();

    data_matrix_ = Eigen::MatrixXd(sample_size_, num_variants_);

    means_.reserve(num_variants_);
    variances_.reserve(num_variants_);
    monomorphic_indices_.reserve(num_variants_ / 100);
}

auto GenotypeLoader::create(
    const std::filesystem::path& bed_path,
    std::shared_ptr<SampleManager> sample_manager)
    -> std::expected<GenotypeLoader, Error>
{
    auto logger = logging::get();
    auto bed_pipe = BedPipe::create(bed_path, std::move(sample_manager));
    if (!bed_pipe)
    {
        return std::unexpected(bed_pipe.error());
    }

    return GenotypeLoader{std::move(*bed_pipe)};
}

auto GenotypeLoader::finalize() -> std::expected<GenotypeMatrix, Error>
{
    // Convert statistics to Eigen vectors
    Eigen::VectorXd mean_vec
        = Eigen::Map<Eigen::VectorXd>(means_.data(), means_.size());
    Eigen::VectorXd variance_vec
        = Eigen::Map<Eigen::VectorXd>(variances_.data(), variances_.size());

    // Convert monomorphic indices to set
    std::unordered_set<int64_t> mono_set(
        monomorphic_indices_.begin(), monomorphic_indices_.end());

    // Create final genotype matrix
    auto matrix = GenotypeMatrix::create(
        std::move(data_matrix_),
        std::move(mono_set),
        std::move(mean_vec),
        std::move(variance_vec));

    return matrix;
}
}  // namespace gelex
