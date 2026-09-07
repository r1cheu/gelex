// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/detail/atomic_output_stream.h"

#include <filesystem>
#include <fmt/format.h>
#include <ios>
#include <string_view>
#include <system_error>
#include <utility>

#include "gelex/exception.h"

namespace gelex::detail
{

namespace
{

auto tmp_path_for(const std::filesystem::path& final_path)
    -> std::filesystem::path
{
    std::filesystem::path tmp = final_path;
    tmp += ".tmp";
    return tmp;
}

}  // namespace

AtomicOutputStream::AtomicOutputStream(
    std::filesystem::path path,
    std::ios::openmode mode)
    : path_(std::move(path)), tmp_path_(tmp_path_for(path_))
{
    if (std::filesystem::is_directory(path_))
    {
        throw GelexException(
            fmt::format(
                "{}: is a directory, not a regular file", path_.string()));
    }

    file_.open(tmp_path_, mode | std::ios::out);
    if (!file_.is_open())
    {
        throw GelexException(
            fmt::format("{}: failed to open file", tmp_path_.string()));
    }
}

AtomicOutputStream::~AtomicOutputStream() noexcept
{
    if (committed_)
    {
        return;
    }

    try
    {
        if (file_.is_open())
        {
            file_.close();
        }
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
    }

    std::error_code ec;
    std::filesystem::remove(tmp_path_, ec);
}

auto AtomicOutputStream::write(const char* data, std::streamsize size) -> void
{
    file_.write(data, size);
    if (!file_)
    {
        throw GelexException(
            fmt::format("{}: failed to write", path_.string()));
    }
}

auto AtomicOutputStream::write(std::string_view text) -> void
{
    write(text.data(), static_cast<std::streamsize>(text.size()));
}

auto AtomicOutputStream::seek(std::streamoff offset) -> void
{
    file_.seekp(offset);
    if (!file_)
    {
        throw GelexException(fmt::format("{}: failed to seek", path_.string()));
    }
}

auto AtomicOutputStream::commit() -> void
{
    if (committed_)
    {
        return;
    }

    file_.close();
    if (!file_)
    {
        throw GelexException(
            fmt::format("{}: failed to close file", tmp_path_.string()));
    }

    std::error_code ec;
    std::filesystem::rename(tmp_path_, path_, ec);
    if (ec)
    {
        throw GelexException(
            fmt::format(
                "{}: failed to rename from \"{}\": {}",
                path_.string(),
                tmp_path_.string(),
                ec.message()));
    }
    committed_ = true;
}

}  // namespace gelex::detail
