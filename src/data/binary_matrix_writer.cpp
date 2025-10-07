#include "binary_matrix_writer.h"

#include <Eigen/Core>
#include <expected>
#include <filesystem>
#include <fstream>
#include <ios>

#include "data/loader.h"
#include "gelex/error.h"

namespace gelex::detail
{

using Eigen::Index;

BinaryMatrixWriter::BinaryMatrixWriter(std::ofstream&& file, std::string&& path)
    : file_(std::move(file)), path_(std::move(path))
{
}

auto BinaryMatrixWriter::create(const std::filesystem::path& file_path)
    -> std::expected<BinaryMatrixWriter, Error>
{
    auto file
        = open_file<std::ofstream>(file_path, std::ios::binary | std::ios::out);

    if (!file)
    {
        return std::unexpected(file.error());
    }
    return BinaryMatrixWriter(std::move(*file), std::move(file_path));
}

auto BinaryMatrixWriter::write(const Eigen::MatrixXd& matrix)
    -> std::expected<void, Error>
{
    if (matrix.size() == 0)
    {
        return {};  // Nothing to write
    }

    // Write matrix data in column-major order (Eigen's default)
    file_.write(
        reinterpret_cast<const char*>(matrix.data()),
        static_cast<std::streamsize>(matrix.size() * sizeof(double)));

    file_.flush();

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to write matrix data to binary file"},
            path_));
    }

    return {};
}

}  // namespace gelex::detail
