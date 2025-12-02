#include "binary_matrix_writer.h"

#include "gelex/exception.h"
#include "parser.h"

namespace gelex::detail
{

BinaryMatrixWriter::BinaryMatrixWriter(const std::filesystem::path& file_path)
    : path_(file_path), io_buffer_(kDefaultBufferSize)
{
    file_ = detail::open_file<std::ofstream>(
        path_, std::ios::binary | std::ios::trunc, io_buffer_);
}

void BinaryMatrixWriter::write(const Eigen::Ref<const Eigen::MatrixXd>& matrix)
{
    if (!file_.good())
    {
        throw FileWriteException(
            std::format(
                "{}:Binary matrix writer stream is in bad state",
                path_.string()));
    }

    if (matrix.size() == 0)
    {
        return;
    }

    file_.write(
        reinterpret_cast<const char*>(matrix.data()),
        static_cast<std::streamsize>(matrix.size() * sizeof(double)));

    if (!file_.good())
    {
        throw FileWriteException(
            std::format(
                "{}:Failed to write matrix data to binary file",
                path_.string()));
    }
}

}  // namespace gelex::detail
