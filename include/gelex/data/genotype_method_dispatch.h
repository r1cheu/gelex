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

#ifndef GELEX_DATA_GENOTYPE_METHOD_DISPATCH_H_
#define GELEX_DATA_GENOTYPE_METHOD_DISPATCH_H_

#include <array>
#include <cctype>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/data/genotype_processor.h"
#include "gelex/exception.h"

namespace gelex
{

enum class GenotypeProcessMethod : uint8_t
{
    Standardize,
    Center,
    OrthStandardize,
    OrthCenter,
    StandardizeSample,
    CenterSample,
    OrthStandardizeSample,
    OrthCenterSample
};

struct GenotypeMethodEntry
{
    std::string_view name;
    GenotypeProcessMethod method;
};

constexpr auto genotype_method_entries = std::array<GenotypeMethodEntry, 8>{
    GenotypeMethodEntry{"standardize", GenotypeProcessMethod::Standardize},
    GenotypeMethodEntry{"center", GenotypeProcessMethod::Center},
    GenotypeMethodEntry{
        "orth-standardize",
        GenotypeProcessMethod::OrthStandardize},
    GenotypeMethodEntry{"orth-center", GenotypeProcessMethod::OrthCenter},
    GenotypeMethodEntry{
        "standardize-sample",
        GenotypeProcessMethod::StandardizeSample},
    GenotypeMethodEntry{"center-sample", GenotypeProcessMethod::CenterSample},
    GenotypeMethodEntry{
        "orth-standardize-sample",
        GenotypeProcessMethod::OrthStandardizeSample},
    GenotypeMethodEntry{
        "orth-center-sample",
        GenotypeProcessMethod::OrthCenterSample}};

constexpr std::string_view genotype_method_hint
    = "standardize, center, orth-standardize, orth-center and with -sample "
      "suffix e.g. standardize-sample";

inline auto to_ascii_lower(std::string_view input) -> std::string
{
    std::string value(input);
    for (auto& ch : value)
    {
        ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    }
    return value;
}

inline auto find_genotype_method_entry_by_name(std::string_view method_name)
    -> const GenotypeMethodEntry*
{
    for (const auto& entry : genotype_method_entries)
    {
        if (entry.name == method_name)
        {
            return &entry;
        }
    }
    return nullptr;
}

inline auto find_genotype_method_entry_by_method(GenotypeProcessMethod method)
    -> const GenotypeMethodEntry*
{
    for (const auto& entry : genotype_method_entries)
    {
        if (entry.method == method)
        {
            return &entry;
        }
    }
    return nullptr;
}

inline auto is_center_family_method(GenotypeProcessMethod method) -> bool
{
    switch (method)
    {
        case GenotypeProcessMethod::Center:
        case GenotypeProcessMethod::OrthCenter:
        case GenotypeProcessMethod::CenterSample:
        case GenotypeProcessMethod::OrthCenterSample:
            return true;
        case GenotypeProcessMethod::Standardize:
        case GenotypeProcessMethod::OrthStandardize:
        case GenotypeProcessMethod::StandardizeSample:
        case GenotypeProcessMethod::OrthStandardizeSample:
            return false;
    }
    throw InvalidInputException("Invalid genotype process method.");
}

inline auto parse_genotype_process_method(std::string_view method)
    -> GenotypeProcessMethod
{
    auto normalized = to_ascii_lower(method);
    const auto* entry = find_genotype_method_entry_by_name(normalized);
    if (entry != nullptr)
    {
        return entry->method;
    }

    throw InvalidInputException(
        "Unknown genotype process method: " + std::string(method)
        + ". Valid: " + std::string(genotype_method_hint));
}

inline auto genotype_process_method_name(GenotypeProcessMethod method)
    -> std::string_view
{
    const auto* entry = find_genotype_method_entry_by_method(method);
    if (entry != nullptr)
    {
        return entry->name;
    }
    throw InvalidInputException(
        "Invalid genotype process method. Valid: "
        + std::string(genotype_method_hint));
}

template <typename Visitor>
inline auto visit_genotype_method(
    GenotypeProcessMethod method,
    Visitor&& visitor) -> decltype(auto)
{
    switch (method)
    {
        case GenotypeProcessMethod::Standardize:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::StandardizedHWE>();
        case GenotypeProcessMethod::Center:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::CenteredHWE>();
        case GenotypeProcessMethod::OrthStandardize:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::OrthStandardizedHWE>();
        case GenotypeProcessMethod::OrthCenter:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::OrthCenteredHWE>();
        case GenotypeProcessMethod::StandardizeSample:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::Standardized>();
        case GenotypeProcessMethod::CenterSample:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::Centered>();
        case GenotypeProcessMethod::OrthStandardizeSample:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::OrthStandardized>();
        case GenotypeProcessMethod::OrthCenterSample:
            return std::forward<Visitor>(visitor)
                .template operator()<grm::OrthCentered>();
    }
    throw InvalidInputException("Invalid genotype process method.");
}

template <typename Visitor>
inline auto visit_assoc_method(GenotypeProcessMethod method, Visitor&& visitor)
    -> decltype(auto)
{
    if (!is_center_family_method(method))
    {
        throw InvalidInputException(
            "assoc --geno-method supports only center-family methods: center, "
            "orth-center, center-sample, orth-center-sample");
    }
    return visit_genotype_method(method, std::forward<Visitor>(visitor));
}

}  // namespace gelex

#endif  // GELEX_DATA_GENOTYPE_METHOD_DISPATCH_H_
