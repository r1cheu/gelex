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

AtomicOutputStream::AtomicOutputStream(std::filesystem::path path)
    : path_(std::move(path)), tmp_path_(tmp_path_for(path_))
{
    if (path_.empty())
    {
        throw GelexException("output path is empty");
    }
    if (std::filesystem::is_directory(path_))
    {
        throw GelexException(
            fmt::format(
                "{}: is a directory, not a regular file", path_.string()));
    }

    file_.open(
        tmp_path_, std::ios::out | std::ios::binary | std::ios::noreplace);

    if (!file_.is_open())
    {
        throw GelexException(
            fmt::format(
                "{}: failed to create temporary file", tmp_path_.string()));
    }
}

AtomicOutputStream::~AtomicOutputStream() noexcept
{
    discard();
}

auto AtomicOutputStream::write(const char* data, std::streamsize size) -> void
{
    file_.write(data, size);
    if (!file_)
    {
        discard();
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
        discard();
        throw GelexException(fmt::format("{}: failed to seek", path_.string()));
    }
}

auto AtomicOutputStream::commit() -> void
{
    if (!file_.is_open())
    {
        throw GelexException(
            fmt::format("{}: file is not open", path_.string()));
    }
    std::error_code ec;
    file_.close();
    if (!file_)
    {
        std::filesystem::remove(tmp_path_, ec);
        throw GelexException(
            fmt::format("{}: failed to close file", path_.string()));
    }
    std::filesystem::rename(tmp_path_, path_, ec);
    if (ec)
    {
        const auto reason = ec.message();
        std::filesystem::remove(tmp_path_, ec);
        throw GelexException(
            fmt::format(
                "{}: failed to rename temporary file: {}",
                path_.string(),
                reason));
    }
}

auto AtomicOutputStream::discard() noexcept -> void
{
    if (!file_.is_open())
    {
        return;
    }
    file_.close();
    std::error_code ec;
    std::filesystem::remove(tmp_path_, ec);
}

}  // namespace gelex::detail
