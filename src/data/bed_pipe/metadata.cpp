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
#include <Eigen/Core>
#include <string>

#include "gelex/data/reader.h"
#include "gelex/io/detail/parser.h"

namespace gelex::genotype::detail
{

auto load_bed_metadata(const std::string& bfile_prefix) -> BedMetadata
{
    BedMetadata metadata;
    metadata.bfile_prefix = bfile_prefix;
    metadata.num_raw_snps = static_cast<Eigen::Index>(
        io::detail::count_total_lines(bfile_prefix + ".bim"));

    metadata.raw_ids = read_fam(bfile_prefix + ".fam").index().take_keys();
    metadata.num_raw_samples
        = static_cast<Eigen::Index>(metadata.raw_ids.size());
    metadata.bytes_per_variant = (metadata.num_raw_samples + 3) / 4;

    return metadata;
}

}  // namespace gelex::genotype::detail
