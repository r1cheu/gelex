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

#ifndef GELEX_INTERNAL_DATA_LOADER_SNP_EFFECT_LOADER_H_
#define GELEX_INTERNAL_DATA_LOADER_SNP_EFFECT_LOADER_H_

#include <algorithm>
#include <filesystem>
#include <span>
#include <string_view>
#include <utility>

#include <Eigen/Core>

#include "gelex/types/snp_info.h"

namespace gelex::detail
{

struct ColumnIndices
{
    int chrom = -1;
    int id = -1;
    int pos = -1;
    int a1 = -1;
    int a2 = -1;
    int a1frq = -1;
    int add = -1;
    int dom = -1;

    [[nodiscard]] bool has_required_columns() const
    {
        return chrom != -1 && id != -1 && pos != -1 && a1 != -1 && a2 != -1
               && a1frq != -1 && add != -1;
    }

    [[nodiscard]] int max_required_index() const
    {
        int m = std::max({chrom, id, pos, a1, a2, a1frq, add});
        if (dom != -1)
        {
            m = std::max(m, dom);
        }
        return m;
    }
};

class SnpEffectLoader
{
   public:
    explicit SnpEffectLoader(const std::filesystem::path& snp_effect_path);

    const SnpEffects& effects() const { return snp_effects_; }
    SnpEffects&& take_effects() && { return std::move(snp_effects_); }
    bool has_dom_effects() const { return has_dom_; }

   private:
    static ColumnIndices assign_column_indices(
        std::span<const std::string_view> header_columns);

    void load(const std::filesystem::path& snp_effect_path);
    static void parse_header(std::string_view line, ColumnIndices& indices);
    void parse_line(
        std::string_view line,
        int line_number,
        const ColumnIndices& indices);

    SnpEffects snp_effects_;
    bool has_dom_ = false;
};

bool check_dom_effect_column(const std::filesystem::path& snp_effect_path);

}  // namespace gelex::detail

#endif  // GELEX_INTERNAL_DATA_LOADER_SNP_EFFECT_LOADER_H_
