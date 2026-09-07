// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_IO_MAPPED_FILE_H_
#define GELEX_IO_MAPPED_FILE_H_

#include <cstddef>
#include <memory>
#include <string>
#include <system_error>

namespace gelex
{

class MappedFile
{
   public:
    MappedFile() noexcept;
    ~MappedFile();
    MappedFile(MappedFile&&) noexcept;
    auto operator=(MappedFile&&) noexcept -> MappedFile&;
    MappedFile(const MappedFile&) = delete;
    auto operator=(const MappedFile&) -> MappedFile& = delete;

    auto map(const std::string& path, std::error_code& ec) -> void;

    [[nodiscard]] auto data() const noexcept -> const std::byte*;
    [[nodiscard]] auto size() const noexcept -> std::size_t;

   private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace gelex

#endif  // GELEX_IO_MAPPED_FILE_H_
