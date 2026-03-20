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

#ifndef GELEX_TYPES_GENOTYPE_PROCESS_METHOD_H_
#define GELEX_TYPES_GENOTYPE_PROCESS_METHOD_H_

#include <cstdint>

#include <fmt/base.h>

#include "gelex/exception.h"

namespace gelex
{

enum class GenotypeProcessMethod : uint8_t
{
    StandardizeHWE = 1,
    CenterHWE,
    OrthStandardizeHWE,
    OrthCenterHWE,
    Standardize,
    Center,
    OrthStandardize,
    OrthCenter
};

inline auto is_center_family_method(GenotypeProcessMethod method) -> bool
{
    switch (method)
    {
        case GenotypeProcessMethod::CenterHWE:
        case GenotypeProcessMethod::OrthCenterHWE:
        case GenotypeProcessMethod::Center:
        case GenotypeProcessMethod::OrthCenter:
            return true;
        case GenotypeProcessMethod::StandardizeHWE:
        case GenotypeProcessMethod::OrthStandardizeHWE:
        case GenotypeProcessMethod::Standardize:
        case GenotypeProcessMethod::OrthStandardize:
            return false;
    }
    throw InvalidInputException("Invalid genotype process method.");
}

inline auto is_orthogonal_method(GenotypeProcessMethod method) -> bool
{
    switch (method)
    {
        case GenotypeProcessMethod::OrthStandardizeHWE:
        case GenotypeProcessMethod::OrthCenterHWE:
        case GenotypeProcessMethod::OrthStandardize:
        case GenotypeProcessMethod::OrthCenter:
            return true;
        case GenotypeProcessMethod::StandardizeHWE:
        case GenotypeProcessMethod::CenterHWE:
        case GenotypeProcessMethod::Standardize:
        case GenotypeProcessMethod::Center:
            return false;
    }
    throw InvalidInputException("Invalid genotype process method.");
}

enum class ModelType : uint8_t
{
    A,
    D,
    AD
};

struct LocusStatistic
{
    double mean{0};
    double stddev{0};
    double maf{0};
    bool is_monomorphic{false};
};

}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::GenotypeProcessMethod> : formatter<string_view>
{
    auto format(gelex::GenotypeProcessMethod t, format_context& ctx) const
        -> format_context::iterator
    {
        string_view name = "unknown";
        using gelex::GenotypeProcessMethod;
        switch (t)
        {
            case GenotypeProcessMethod::StandardizeHWE:
                name = "StandardizeHWE";
                break;
            case GenotypeProcessMethod::CenterHWE:
                name = "CenterHWE";
                break;
            case GenotypeProcessMethod::OrthStandardizeHWE:
                name = "OrthStandardizeHWE";
                break;
            case GenotypeProcessMethod::OrthCenterHWE:
                name = "OrthCenterHWE";
                break;
            case GenotypeProcessMethod::Standardize:
                name = "Standardize";
                break;
            case GenotypeProcessMethod::Center:
                name = "Center";
                break;
            case GenotypeProcessMethod::OrthStandardize:
                name = "OrthStandardize";
                break;
            case GenotypeProcessMethod::OrthCenter:
                name = "OrthCenter";
                break;
            default:
                name = "Unknown Genotype Process Method";
                break;
        }

        return formatter<string_view>::format(name, ctx);
    }
};
}  // namespace fmt

#endif  // GELEX_TYPES_GENOTYPE_PROCESS_METHOD_H_
