#ifndef GELEX_DATA_GRM_BIN_WRITER_H
#define GELEX_DATA_GRM_BIN_WRITER_H

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

namespace gelex::detail
{

class GrmBinWriter
{
   public:
    static constexpr size_t kDefaultBufferSize = static_cast<size_t>(64 * 1024);

    explicit GrmBinWriter(const std::filesystem::path& file_path);

    GrmBinWriter(const GrmBinWriter&) = delete;
    GrmBinWriter(GrmBinWriter&&) noexcept = default;
    GrmBinWriter& operator=(const GrmBinWriter&) = delete;
    GrmBinWriter& operator=(GrmBinWriter&&) noexcept = default;
    ~GrmBinWriter() = default;

    // Write GRM matrix (unnormalized)
    // Format: [float32 lower triangle]
    // Order: (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
    auto write(const Eigen::Ref<const Eigen::MatrixXd>& grm) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::vector<char> io_buffer_;
    std::ofstream file_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_GRM_BIN_WRITER_H
