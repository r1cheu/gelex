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
