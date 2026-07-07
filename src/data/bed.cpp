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

#include "gelex/data/bed.h"

#include <filesystem>
#include <string>
#include <utility>

#include "gelex/data/detail/bed_source.h"
#include "gelex/data/reader.h"

namespace gelex
{

auto open_bed(const std::string& bfile_prefix) -> Bed
{
    auto source_index
        = read_fam(std::filesystem::path{bfile_prefix + ".fam"}).index();
    auto bim = read_bim(std::filesystem::path{bfile_prefix + ".bim"});
    auto bed_source = detail::open_bed_source(bfile_prefix);

    return Bed{std::move(bed_source), std::move(source_index), std::move(bim)};
}

}  // namespace gelex
