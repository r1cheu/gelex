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

#include "gelex/io/detail/atomic_ofstream.h"

#include <exception>
#include <filesystem>
#include <ios>
#include <span>
#include <system_error>
#include <utility>

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/infra/logger.h"

namespace gelex::io::detail
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

AtomicOfstream::AtomicOfstream(
    std::filesystem::path final_path,
    std::ios::openmode mode,
    std::span<char> custom_buffer)
    : final_path_(std::move(final_path)), tmp_path_(tmp_path_for(final_path_))
{
    if (std::filesystem::is_directory(final_path_))
    {
        throw GelexException(
            fmt::format(
                "{}: is a directory, not a regular file",
                final_path_.string()));
    }

    if (!custom_buffer.empty())
    {
        rdbuf()->pubsetbuf(
            custom_buffer.data(),
            static_cast<std::streamsize>(custom_buffer.size()));
    }

    open(tmp_path_, mode | std::ios::out);
    if (!is_open())
    {
        throw GelexException(
            fmt::format("{}: failed to open file", tmp_path_.string()));
    }
}

AtomicOfstream::~AtomicOfstream() noexcept
{
    if (committed_)
    {
        return;
    }

    const bool keep = good() && std::uncaught_exceptions() == 0;

    try
    {
        close();
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
    }

    std::error_code ec;
    if (keep)
    {
        std::filesystem::rename(tmp_path_, final_path_, ec);
        if (!ec)
        {
            return;
        }
        if (auto logger = gelex::logging::get())
        {
            logger->error(
                "{}: failed to rename from \"{}\", discarding output: {}",
                final_path_.string(),
                tmp_path_.string(),
                ec.message());
        }
    }
    std::filesystem::remove(tmp_path_, ec);
}

auto AtomicOfstream::commit() -> void
{
    if (committed_)
    {
        return;
    }
    if (!good())
    {
        throw GelexException(
            fmt::format("{}: stream is not writable", final_path_.string()));
    }

    close();
    if (!good())
    {
        throw GelexException(
            fmt::format("{}: failed to close file", tmp_path_.string()));
    }

    std::error_code ec;
    std::filesystem::rename(tmp_path_, final_path_, ec);
    if (ec)
    {
        throw GelexException(
            fmt::format(
                "{}: failed to rename from \"{}\": {}",
                final_path_.string(),
                tmp_path_.string(),
                ec.message()));
    }
    committed_ = true;
}

}  // namespace gelex::io::detail
