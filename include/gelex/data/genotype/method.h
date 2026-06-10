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

#ifndef GELEX_DATA_GENOTYPE_METHOD_H_
#define GELEX_DATA_GENOTYPE_METHOD_H_

#include <cstdint>
#include <utility>

#include <fmt/base.h>
#include <fmt/format.h>

#include "gelex/exception.h"

namespace gelex
{

enum class GenotypeMethod : uint8_t
{
    StandardizeHWE = 0,
    CenterHWE = 1,
    Standardize = 2,
    Center = 3,
    OrthStandardizeHWE = 4,
    OrthCenterHWE = 5,
    OrthStandardize = 6,
    OrthCenter = 7,
    NOIAStandardize = 8,
    NOIACenter = 9
};

constexpr auto is_center(GenotypeMethod method) -> bool
{
    switch (method)
    {
        case GenotypeMethod::CenterHWE:
        case GenotypeMethod::Center:
        case GenotypeMethod::OrthCenterHWE:
        case GenotypeMethod::OrthCenter:
        case GenotypeMethod::NOIACenter:
            return true;
        default:
            return false;
    }
}

constexpr auto is_orthogonal(GenotypeMethod method) -> bool
{
    switch (method)
    {
        case GenotypeMethod::OrthStandardizeHWE:
        case GenotypeMethod::OrthCenterHWE:
        case GenotypeMethod::OrthStandardize:
        case GenotypeMethod::OrthCenter:
            return true;
        default:
            return false;
    }
}

constexpr auto is_noia(GenotypeMethod method) -> bool
{
    switch (method)
    {
        case GenotypeMethod::NOIAStandardize:
        case GenotypeMethod::NOIACenter:
            return true;
        default:
            return false;
    }
}

constexpr auto is_hwe(GenotypeMethod method) -> bool
{
    switch (method)
    {
        case GenotypeMethod::StandardizeHWE:
        case GenotypeMethod::CenterHWE:
        case GenotypeMethod::OrthStandardizeHWE:
        case GenotypeMethod::OrthCenterHWE:
            return true;
        default:
            return false;
    }
}

inline auto genotype_method_from_byte(uint8_t b) -> GenotypeMethod
{
    switch (b)
    {
        case std::to_underlying(GenotypeMethod::StandardizeHWE):
            return GenotypeMethod::StandardizeHWE;
        case std::to_underlying(GenotypeMethod::CenterHWE):
            return GenotypeMethod::CenterHWE;
        case std::to_underlying(GenotypeMethod::Standardize):
            return GenotypeMethod::Standardize;
        case std::to_underlying(GenotypeMethod::Center):
            return GenotypeMethod::Center;
        case std::to_underlying(GenotypeMethod::OrthStandardizeHWE):
            return GenotypeMethod::OrthStandardizeHWE;
        case std::to_underlying(GenotypeMethod::OrthCenterHWE):
            return GenotypeMethod::OrthCenterHWE;
        case std::to_underlying(GenotypeMethod::OrthStandardize):
            return GenotypeMethod::OrthStandardize;
        case std::to_underlying(GenotypeMethod::OrthCenter):
            return GenotypeMethod::OrthCenter;
        case std::to_underlying(GenotypeMethod::NOIAStandardize):
            return GenotypeMethod::NOIAStandardize;
        case std::to_underlying(GenotypeMethod::NOIACenter):
            return GenotypeMethod::NOIACenter;
        default:
            throw GelexException(
                fmt::format("Invalid genotype method byte: {}", b));
    }
}

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
struct formatter<gelex::GenotypeMethod> : formatter<string_view>
{
    static auto format(gelex::GenotypeMethod method, format_context& ctx)
        -> format_context::iterator
    {
        string_view name = "Unknown";
        switch (method)
        {
            case gelex::GenotypeMethod::StandardizeHWE:
                name = "StandardizeHWE";
                break;
            case gelex::GenotypeMethod::CenterHWE:
                name = "CenterHWE";
                break;
            case gelex::GenotypeMethod::Standardize:
                name = "Standardize";
                break;
            case gelex::GenotypeMethod::Center:
                name = "Center";
                break;
            case gelex::GenotypeMethod::OrthStandardizeHWE:
                name = "OrthStandardizeHWE";
                break;
            case gelex::GenotypeMethod::OrthCenterHWE:
                name = "OrthCenterHWE";
                break;
            case gelex::GenotypeMethod::OrthStandardize:
                name = "OrthStandardize";
                break;
            case gelex::GenotypeMethod::OrthCenter:
                name = "OrthCenter";
                break;
            case gelex::GenotypeMethod::NOIAStandardize:
                name = "NOIAStandardize";
                break;
            case gelex::GenotypeMethod::NOIACenter:
                name = "NOIACenter";
                break;
        }
        return fmt::format_to(ctx.out(), "{}", name);
    }
};
}  // namespace fmt

#endif  // GELEX_DATA_GENOTYPE_METHOD_H_
