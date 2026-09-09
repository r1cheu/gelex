// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DENSE_READER_H_
#define GELEX_IO_DENSE_READER_H_

#include <Eigen/Core>
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

class DenseReader
{
   public:
    explicit DenseReader(std::string_view file_path);

    DenseReader(const DenseReader&) = delete;
    DenseReader(DenseReader&&) noexcept = default;
    auto operator=(const DenseReader&) -> DenseReader& = delete;
    auto operator=(DenseReader&&) noexcept -> DenseReader& = default;
    ~DenseReader() = default;

    auto contains(std::string_view identifier) const -> bool;
    auto info(std::string_view identifier) const -> const MatrixHeader&;
    auto payloads() const -> std::vector<MatrixHeader>;

    // The map aliases the memory-mapped file and must not outlive the reader.
    template <detail::SupportedDtype T>
    auto to_map(std::string_view identifier) const
        -> Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>;

    template <detail::SupportedDtype T>
    auto to_mat(std::string_view identifier) const -> Eigen::MatrixX<T>;

    auto size() const noexcept -> std::size_t;

   private:
    auto parse_footer_and_index() -> void;

    auto find_entry(std::string_view identifier) const
        -> const detail::MatrixEntry&;

    auto payload_bytes(const detail::MatrixEntry& entry) const
        -> std::span<const std::byte>;

    auto validate_payload_ranges(
        std::vector<std::pair<std::uint64_t, std::uint64_t>> ranges) const
        -> void;

    std::filesystem::path path_;
    MappedFile mmap_;
    std::vector<detail::MatrixEntry> payloads_;
    std::unordered_map<std::string_view, std::size_t> index_;
};

template <detail::SupportedDtype T>
auto DenseReader::to_map(std::string_view identifier) const
    -> Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>
{
    const auto& entry = find_entry(identifier);
    const auto& header = entry.header;
    if (header.type != detail::binary_type_for<T>)
    {
        throw GelexException(
            fmt::format(
                "{}: dtype mismatch for payload \"{}\": stored={}, "
                "requested={}",
                path_.string(),
                identifier,
                std::to_underlying(header.type),
                std::to_underlying(detail::binary_type_for<T>)));
    }
    if (!std::in_range<Eigen::Index>(header.shape[0])
        || !std::in_range<Eigen::Index>(header.shape[1]))
    {
        throw GelexException(
            fmt::format(
                "{}: payload \"{}\" shape exceeds Eigen::Index",
                path_.string(),
                identifier));
    }

    const auto bytes = payload_bytes(entry);
    const auto* data = reinterpret_cast<const T*>(bytes.data());
    return Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>(
        data,
        static_cast<Eigen::Index>(header.shape[0]),
        static_cast<Eigen::Index>(header.shape[1]));
}

template <detail::SupportedDtype T>
auto DenseReader::to_mat(std::string_view identifier) const -> Eigen::MatrixX<T>
{
    return Eigen::MatrixX<T>{to_map<T>(identifier)};
}

}  // namespace gelex

#endif  // GELEX_IO_DENSE_READER_H_
