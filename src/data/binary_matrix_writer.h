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
    static auto create(const std::filesystem::path& file_path)
        -> std::expected<BinaryMatrixWriter, Error>;

    BinaryMatrixWriter(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter(BinaryMatrixWriter&&) noexcept = default;
    BinaryMatrixWriter& operator=(const BinaryMatrixWriter&) = delete;
    BinaryMatrixWriter& operator=(BinaryMatrixWriter&&) noexcept = default;
    ~BinaryMatrixWriter();

    auto write(const Eigen::MatrixXd& matrix) -> std::expected<void, Error>;

   private:
    explicit BinaryMatrixWriter(std::ofstream&& file, std::string&& path);
    std::ofstream file_;
    std::filesystem::path path_;
};

}  // namespace gelex::detail
