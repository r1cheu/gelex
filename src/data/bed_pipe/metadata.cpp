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

#include "metadata.h"

#include "gelex/internal/data/loader/fam_loader.h"
#include "gelex/io/parser.h"

namespace gelex::detail
{

auto load_bed_metadata(const std::filesystem::path& bed_prefix) -> BedMetadata
{
    auto bed_path = bed_prefix;
    bed_path.replace_extension(".bed");
    auto bim_path = bed_prefix;
    bim_path.replace_extension(".bim");
    auto fam_path = bed_prefix;
    fam_path.replace_extension(".fam");

    BedMetadata metadata;
    metadata.bed_path = bed_path;
    metadata.num_raw_snps
        = static_cast<Eigen::Index>(count_total_lines(bim_path));

    FamLoader fam_loader(fam_path);
    metadata.raw_ids = fam_loader.ids();
    metadata.num_raw_samples
        = static_cast<Eigen::Index>(metadata.raw_ids.size());
    metadata.bytes_per_variant = (metadata.num_raw_samples + 3) / 4;

    return metadata;
}

}  // namespace gelex::detail
