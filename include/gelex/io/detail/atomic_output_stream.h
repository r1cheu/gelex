// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
#define GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_

#include <concepts>
#include <filesystem>
#include <fstream>
#include <string_view>

namespace gelex::detail
{

class AtomicOutputStream
{
   public:
    explicit AtomicOutputStream(std::filesystem::path path);

    AtomicOutputStream(const AtomicOutputStream&) = delete;
    // The moved-from stream is closed and owns nothing.
    AtomicOutputStream(AtomicOutputStream&&) noexcept = default;
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
    auto discard() noexcept -> void;

    std::filesystem::path path_;
    std::filesystem::path tmp_path_;
    std::ofstream file_;
};

template <std::unsigned_integral T>
auto write_integer(AtomicOutputStream& file, T value) -> void
{
    file.write(
        reinterpret_cast<const char*>(&value),
        static_cast<std::streamsize>(sizeof(value)));
}

}  // namespace gelex::detail

#endif  // GELEX_IO_DETAIL_ATOMIC_OUTPUT_STREAM_H_
