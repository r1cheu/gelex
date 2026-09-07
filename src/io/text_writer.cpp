// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/detail/text_writer.h"

#include <exception>
#include <filesystem>
#include <fmt/format.h>
#include <initializer_list>
#include <ios>
#include <string_view>

#include "gelex/infra/log.h"

namespace gelex::detail
{

TextWriter::TextWriter(const std::filesystem::path& path)
    : ofs_(path, std::ios::out)
{
}

TextWriter::~TextWriter() noexcept
{
    if (std::uncaught_exceptions() > 0)
    {
        return;
    }
    try
    {
        ofs_.commit();
    }
    catch (const std::exception& exception)
    {
        try
        {
            error(
                fmt::format(
                    "{}: failed to commit, discarding output: {}",
                    ofs_.path().string(),
                    exception.what()));
        }
        catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
        {
        }
    }
    catch (...)  // NOLINT(bugprone-empty-catch): dtor must be noexcept
    {
    }
}

auto TextWriter::write_header(std::initializer_list<std::string_view> columns)
    -> void
{
    bool first = true;
    for (const auto col : columns)
    {
        if (!first)
        {
            ofs_.write("\t");
        }
        ofs_.write(col);
        first = false;
    }
    ofs_.write("\n");
}

auto TextWriter::write(std::string_view line) -> void
{
    ofs_.write(line);
    ofs_.write("\n");
}

auto TextWriter::path() const noexcept -> const std::filesystem::path&
{
    return ofs_.path();
}

}  // namespace gelex::detail
