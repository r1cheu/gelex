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

#ifndef GELEX_TYPES_CHR_GROUP_H_
#define GELEX_TYPES_CHR_GROUP_H_

#include <string>
#include <utility>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class SnpEffects;

struct ChrGroup
{
    std::string name;
    std::vector<std::pair<Eigen::Index, Eigen::Index>> ranges;
    Eigen::Index total_snps;
};

auto build_chr_groups(bool do_loco, const SnpEffects& snp_effects)
    -> std::vector<ChrGroup>;

}  // namespace gelex

#endif  // GELEX_TYPES_CHR_GROUP_H_
