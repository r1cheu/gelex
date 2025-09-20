#pragma once

#include <expected>
#include <filesystem>
#include <fstream>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex::detail
{

class BinaryMatrixWriter
{
   public:
    explicit BinaryMatrixWriter(std::filesystem::path file_path);

    BinaryMatrixWriter(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter(BinaryMatrixWriter&&) noexcept = default;
    BinaryMatrixWriter& operator=(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter& operator=(BinaryMatrixWriter&&) noexcept = default;
    ~BinaryMatrixWriter();

    auto open() -> std::expected<void, Error>;
    auto close() -> std::expected<void, Error>;

    auto append_matrix(const Eigen::MatrixXd& matrix)
        -> std::expected<void, Error>;

    const std::filesystem::path& file_path() const noexcept
    {
        return file_path_;
    }
    bool is_open() const noexcept { return file_.is_open(); }

   private:
    std::filesystem::path file_path_;
    std::ofstream file_;
};

}  // namespace gelex::detail
