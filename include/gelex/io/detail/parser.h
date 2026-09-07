// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DETAIL_PARSER_H_
#define GELEX_IO_DETAIL_PARSER_H_

#include <concepts>
#include <cstddef>
#include <filesystem>
#include <fmt/format.h>
#include <ios>
#include <span>
#include <system_error>

#include "gelex/exception.h"

namespace gelex::detail
{

template <typename T>
concept FileStream
    = std::derived_from<T, std::ios>
      && requires(std::filesystem::path p, std::ios::openmode m) {
             { T() };
         };

template <FileStream StreamType>
[[nodiscard]] StreamType open_file(
    const std::filesystem::path& path,
    std::ios::openmode mode,
    std::span<char> custom_buffer = {})
{
    if (std::filesystem::is_directory(path))
    {
        throw gelex::GelexException(
            fmt::format(
                "{}: is a directory, not a regular file", path.string()));
    }

    StreamType stream;
    if (!custom_buffer.empty())
    {
        stream.rdbuf()->pubsetbuf(
            custom_buffer.data(),
            static_cast<std::streamsize>(custom_buffer.size()));
    }

    stream.open(path, mode);

    if (!stream.is_open())
    {
        if ((mode & std::ios::in) && !std::filesystem::exists(path))
        {
            throw GelexException(fmt::format("{}: not found", path.string()));
        }
        throw GelexException(
            fmt::format("{}: failed to open file", path.string()));
    }

    if ((mode & std::ios::in) && std::filesystem::is_regular_file(path))
    {
        std::error_code ec;
        if (std::filesystem::file_size(path, ec) == 0 && !ec)
        {
            throw GelexException(fmt::format("{}: is empty", path.string()));
        }
    }

    return stream;
}

size_t count_total_lines(const std::filesystem::path& path);

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_PARSER_H_
