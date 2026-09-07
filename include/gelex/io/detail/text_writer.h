// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DETAIL_TEXT_WRITER_H_
#define GELEX_IO_DETAIL_TEXT_WRITER_H_

#include <filesystem>
#include <initializer_list>
#include <string_view>

#include "gelex/io/detail/atomic_output_stream.h"

namespace gelex::detail
{

class TextWriter
{
   public:
    explicit TextWriter(const std::filesystem::path& path);
    TextWriter(const TextWriter&) = delete;
    TextWriter(TextWriter&&) = delete;
    auto operator=(const TextWriter&) -> TextWriter& = delete;
    auto operator=(TextWriter&&) -> TextWriter& = delete;
    ~TextWriter() noexcept;

    auto write_header(std::initializer_list<std::string_view> columns) -> void;
    auto write(std::string_view line) -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&;

   private:
    AtomicOutputStream ofs_;
};

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_TEXT_WRITER_H_
