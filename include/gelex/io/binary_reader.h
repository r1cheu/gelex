/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_IO_BINARY_READER_H_
#define GELEX_IO_BINARY_READER_H_

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

class BinaryReader
{
   public:
    explicit BinaryReader(std::string_view file_path);

    BinaryReader(const BinaryReader&) = delete;
    BinaryReader(BinaryReader&&) noexcept = default;
    auto operator=(const BinaryReader&) -> BinaryReader& = delete;
    auto operator=(BinaryReader&&) noexcept -> BinaryReader& = default;
    ~BinaryReader() = default;

    [[nodiscard]] auto contains(std::string_view identifier) const -> bool;
    [[nodiscard]] auto info(std::string_view identifier) const& -> const
        PayloadInfo&;
    [[nodiscard]] auto info(std::string_view identifier) const&& -> const
        PayloadInfo& = delete;
    [[nodiscard]] auto payloads() const -> std::vector<PayloadInfo>;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto to_map(std::string_view identifier)
        const& -> Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto to_map(std::string_view identifier)
        const&& -> Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>
        = delete;

    template <detail::SupportedDtype T>
    [[nodiscard]] auto to_mat(std::string_view identifier) const
        -> Eigen::MatrixX<T>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;

   private:
    auto parse_footer_and_directory() -> void;

    [[nodiscard]] auto find_entry(std::string_view identifier) const
        -> const detail::PayloadEntry&;

    [[nodiscard]] auto payload_bytes(const detail::PayloadEntry& entry) const
        -> std::span<const std::byte>;

    auto validate_payload_ranges(
        std::vector<std::pair<std::uint64_t, std::uint64_t>> ranges) const
        -> void;

    std::filesystem::path path_;
    MappedFile mmap_;
    std::vector<detail::PayloadEntry> payloads_;
    std::unordered_map<std::string_view, std::size_t> index_;
};

template <detail::SupportedDtype T>
auto BinaryReader::to_map(std::string_view identifier)
    const& -> Eigen::Map<const Eigen::MatrixX<T>, Eigen::Aligned64>
{
    const auto& entry = find_entry(identifier);
    const auto& descriptor = entry.info.descriptor;
    if (descriptor.type != detail::binary_type_for<T>)
    {
        throw GelexException(
            fmt::format(
                "{}: dtype mismatch for payload \"{}\": stored={}, "
                "requested={}",
                path_.string(),
                identifier,
                std::to_underlying(descriptor.type),
                std::to_underlying(detail::binary_type_for<T>)));
    }
    if (!std::in_range<Eigen::Index>(descriptor.shape[0])
        || !std::in_range<Eigen::Index>(descriptor.shape[1]))
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
        static_cast<Eigen::Index>(descriptor.shape[0]),
        static_cast<Eigen::Index>(descriptor.shape[1]));
}

template <detail::SupportedDtype T>
auto BinaryReader::to_mat(std::string_view identifier) const
    -> Eigen::MatrixX<T>
{
    return Eigen::MatrixX<T>{to_map<T>(identifier)};
}

}  // namespace gelex

#endif  // GELEX_IO_BINARY_READER_H_
