// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_CSC_READER_H_
#define GELEX_IO_CSC_READER_H_

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fmt/format.h>
#include <span>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/detail/binary_wire.h"
#include "gelex/io/mapped_file.h"

namespace gelex
{

// Memory-mapped reader for GELEXSC1 files written by CscWriter.
class CscReader
{
   public:
    template <detail::SupportedDtype T>
    using sparse_matrix_type
        = Eigen::SparseMatrix<T, Eigen::ColMajor, std::int64_t>;
    template <detail::SupportedDtype T>
    using sparse_map_type = Eigen::Map<const sparse_matrix_type<T>>;

    explicit CscReader(std::string_view file_path);

    CscReader(const CscReader&) = delete;
    CscReader(CscReader&&) noexcept = default;
    auto operator=(const CscReader&) -> CscReader& = delete;
    auto operator=(CscReader&&) noexcept -> CscReader& = default;
    ~CscReader() = default;

    auto contains(std::string_view identifier) const -> bool;
    auto info(std::string_view identifier) const -> const MatrixHeader&;
    auto nnz(std::string_view identifier) const -> std::uint64_t;
    auto payloads() const -> std::vector<MatrixHeader>;

    // The map aliases the memory-mapped file and must not outlive the reader.
    template <detail::SupportedDtype T>
    auto to_map(std::string_view identifier) const -> sparse_map_type<T>;

    // Copies the matrix out of the file.
    template <detail::SupportedDtype T>
    auto to_mat(std::string_view identifier) const -> sparse_matrix_type<T>;

    auto size() const noexcept -> std::size_t;

   private:
    auto parse_footer_and_index() -> void;

    auto find_entry(std::string_view identifier) const
        -> const detail::CscEntry&;

    template <typename T>
    auto array_at(std::uint64_t offset) const -> const T*
    {
        const std::span<const std::byte> file{mmap_.data(), mmap_.size()};
        return reinterpret_cast<const T*>(
            file.subspan(static_cast<std::size_t>(offset)).data());
    }

    std::filesystem::path path_;
    MappedFile mmap_;
    std::vector<detail::CscEntry> entries_;
    std::unordered_map<std::string_view, std::size_t> index_;
};

template <detail::SupportedDtype T>
auto CscReader::to_map(std::string_view identifier) const -> sparse_map_type<T>
{
    const auto& entry = find_entry(identifier);
    const auto& header = entry.header;
    if (header.type != detail::binary_type_for<T>)
    {
        throw GelexException(
            fmt::format(
                "{}: dtype mismatch for matrix \"{}\": stored={}, "
                "requested={}",
                path_.string(),
                identifier,
                std::to_underlying(header.type),
                std::to_underlying(detail::binary_type_for<T>)));
    }
    if (!std::in_range<Eigen::Index>(header.shape[0])
        || !std::in_range<Eigen::Index>(header.shape[1])
        || !std::in_range<Eigen::Index>(entry.nnz))
    {
        throw GelexException(
            fmt::format(
                "{}: matrix \"{}\" shape exceeds Eigen::Index",
                path_.string(),
                identifier));
    }
    return sparse_map_type<T>{
        static_cast<Eigen::Index>(header.shape[0]),
        static_cast<Eigen::Index>(header.shape[1]),
        static_cast<Eigen::Index>(entry.nnz),
        array_at<std::int64_t>(entry.offsets[2]),
        array_at<std::int64_t>(entry.offsets[1]),
        array_at<T>(entry.offsets[0])};
}

template <detail::SupportedDtype T>
auto CscReader::to_mat(std::string_view identifier) const
    -> sparse_matrix_type<T>
{
    return sparse_matrix_type<T>{to_map<T>(identifier)};
}

}  // namespace gelex

#endif  // GELEX_IO_CSC_READER_H_
