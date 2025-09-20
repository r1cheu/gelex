#include "binary_matrix_writer.h"

#include <expected>
#include <filesystem>
#include <fstream>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex::detail
{

using Eigen::Index;

BinaryMatrixWriter::BinaryMatrixWriter(std::filesystem::path file_path)
    : file_path_(std::move(file_path))
{
}

BinaryMatrixWriter::~BinaryMatrixWriter()
{
    if (file_.is_open())
    {
        file_.close();
    }
}

auto BinaryMatrixWriter::open() -> std::expected<void, Error>
{
    file_.open(file_path_, std::ios::binary | std::ios::app);
    if (!file_.is_open())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to open binary matrix file"},
            file_path_.string()));
    }

    return {};
}

auto BinaryMatrixWriter::close() -> std::expected<void, Error>
{
    if (file_.is_open())
    {
        file_.close();
    }

    if (file_.fail())
    {
        return std::unexpected(enrich_with_file_info(
            Error{ErrorCode::FileIOError, "Failed to close binary matrix file"},
            file_path_.string()));
    }

    return {};
}

auto BinaryMatrixWriter::append_matrix(const Eigen::MatrixXd& matrix)
    -> std::expected<void, Error>
{
    if (!file_.is_open())
    {
        return std::unexpected(
            Error{ErrorCode::FileIOError, "Binary matrix file is not open"});
    }

    if (matrix.size() == 0)
    {
        return {};  // Nothing to write
    }

    // Write matrix data in column-major order (Eigen's default)
    file_.write(
        reinterpret_cast<const char*>(matrix.data()),
        matrix.size() * sizeof(double));

    if (!file_.good())
    {
        return std::unexpected(enrich_with_file_info(
            Error{
                ErrorCode::FileIOError,
                "Failed to write matrix data to binary file"},
            file_path_.string()));
    }

    return {};
}

}  // namespace gelex::detail
