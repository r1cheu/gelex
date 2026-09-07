// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
#define GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_

#include <filesystem>
#include <fstream>
#include <ios>
#include <string_view>

namespace gelex::detail
{

class AtomicOutputStream
{
   public:
    AtomicOutputStream(std::filesystem::path path, std::ios::openmode mode);

    AtomicOutputStream(const AtomicOutputStream&) = delete;
    AtomicOutputStream(AtomicOutputStream&&) = delete;
    auto operator=(const AtomicOutputStream&) -> AtomicOutputStream& = delete;
    auto operator=(AtomicOutputStream&&) -> AtomicOutputStream& = delete;

    ~AtomicOutputStream() noexcept;

    auto write(const char* data, std::streamsize size) -> void;
    auto write(std::string_view text) -> void;
    auto seek(std::streamoff offset) -> void;
    auto commit() -> void;

    [[nodiscard]] auto path() const noexcept -> const std::filesystem::path&
    {
        return path_;
    }

   private:
    std::filesystem::path path_;
    std::filesystem::path tmp_path_;
    std::ofstream file_;
    bool committed_{false};
};

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
