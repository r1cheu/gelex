#pragma once

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

#include <Eigen/Core>

namespace gelex::detail
{

class BinaryMatrixWriter
{
   public:
    static constexpr size_t kDefaultBufferSize
        = static_cast<const size_t>(64 * 1024);

    explicit BinaryMatrixWriter(const std::filesystem::path& file_path);

    BinaryMatrixWriter(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter(BinaryMatrixWriter&&) noexcept = default;
    BinaryMatrixWriter& operator=(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter& operator=(BinaryMatrixWriter&&) noexcept = default;
    ~BinaryMatrixWriter() = default;

    void write(const Eigen::Ref<const Eigen::MatrixXd>& matrix);

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
