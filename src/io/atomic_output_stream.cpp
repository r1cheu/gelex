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

#include "gelex/io/detail/atomic_output_stream.h"

#include <filesystem>
#include <ios>
#include <string_view>
#include <system_error>
#include <utility>

#include <fmt/format.h>

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
