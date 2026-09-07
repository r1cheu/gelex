// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include "gelex/io/mapped_file.h"

#include <cstddef>
#include <memory>
#include <mio.h>
#include <string>
#include <system_error>

namespace gelex
{

struct MappedFile::Impl
{
    mio::mmap_source mmap;
};

MappedFile::MappedFile() noexcept = default;
MappedFile::~MappedFile() = default;
MappedFile::MappedFile(MappedFile&&) noexcept = default;
auto MappedFile::operator=(MappedFile&&) noexcept -> MappedFile& = default;

auto MappedFile::map(const std::string& path, std::error_code& ec) -> void
{
    if (!impl_)
    {
        impl_ = std::make_unique<Impl>();
    }
    impl_->mmap.map(path, ec);
}

auto MappedFile::data() const noexcept -> const std::byte*
{
    if (!impl_)
    {
        return nullptr;
    }
    return reinterpret_cast<const std::byte*>(impl_->mmap.data());
}

auto MappedFile::size() const noexcept -> std::size_t
{
    return impl_ ? impl_->mmap.size() : 0;
}

}  // namespace gelex
