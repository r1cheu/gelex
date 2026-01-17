#ifndef GELEX_DATA_GENOTYPE_MMAP_H_
#define GELEX_DATA_GENOTYPE_MMAP_H_

#include <filesystem>
#include <fstream>
#include <vector>

#include <mio.h>

#include <Eigen/Core>
#include "gelex/exception.h"

namespace gelex
{

#ifdef USE_AVX512
static constexpr int MAP_OPTIONS = Eigen::Aligned64;
static constexpr size_t ALIGNMENT_BYTES = 64;
#else
static constexpr int MAP_OPTIONS = Eigen::Aligned32;
static constexpr size_t ALIGNMENT_BYTES = 32;
#endif

namespace detail
{
template <typename T>
T read_scalar(std::ifstream& ifs, std::string_view context)
{
    T value;
    ifs.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!ifs)
    {
        throw FileOpenException(
            std::format("Failed to read scalar: {}", context));
    }
    return value;
}

template <typename T>
void read_binary(
    std::ifstream& ifs,
    T* dest,
    size_t count,
    std::string_view context)
{
    if (count == 0)
    {
        return;
    }
    ifs.read(reinterpret_cast<char*>(dest), sizeof(T) * count);
    if (!ifs)
    {
        throw FileOpenException(
            std::format("Failed to read {} data.", context));
    }
}
}  // namespace detail
class GenotypeMap
{
   public:
    using MatrixType = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

    using MapType = Eigen::Map<const Eigen::MatrixXd, MAP_OPTIONS>;

    explicit GenotypeMap(const std::filesystem::path& bin_file);

    GenotypeMap(const GenotypeMap&) = delete;
    GenotypeMap& operator=(const GenotypeMap&) = delete;
    GenotypeMap(GenotypeMap&&) noexcept = default;
    GenotypeMap& operator=(GenotypeMap&&) noexcept = default;
    ~GenotypeMap() = default;

    [[nodiscard]] const MapType& matrix() const noexcept { return mat_; }

    [[nodiscard]] bool is_monomorphic(Eigen::Index snp_index) const noexcept;

    [[nodiscard]] const Eigen::VectorXd& mean() const noexcept { return mean_; }
    [[nodiscard]] const Eigen::VectorXd& stddev() const noexcept
    {
        return stddev_;
    }

    [[nodiscard]] int64_t num_mono() const noexcept
    {
        return static_cast<int64_t>(mono_indices_.size());
    }
    [[nodiscard]] int64_t rows() const noexcept { return rows_; }
    [[nodiscard]] int64_t cols() const noexcept { return cols_; }

   private:
    mio::mmap_source mmap_;

    MapType mat_;

    std::vector<int64_t> mono_indices_;
    Eigen::VectorXd mean_;
    Eigen::VectorXd stddev_;

    int64_t rows_{0};
    int64_t cols_{0};

    void load_metadata(const std::filesystem::path& meta_path);
    static void validate_alignment(const void* ptr);
};

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_MMAP_H_
