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
#include <fmt/format.h>

namespace gelex
{

enum class ScaleMethod : uint8_t
{
    Standardize,
    Center
};

enum class FreqSource : uint8_t
{
    HWE,
    Sample
};

enum class Projection : uint8_t
{
    Direct,
    Orthogonal
};

struct GenotypeProcessMethod
{
    ScaleMethod scale;
    FreqSource freq;
    Projection proj;

    auto operator==(const GenotypeProcessMethod&) const -> bool = default;

    constexpr auto is_center() const -> bool
    {
        return scale == ScaleMethod::Center;
    }
    constexpr auto is_orthogonal() const -> bool
    {
        return proj == Projection::Orthogonal;
    }
    constexpr auto is_hwe() const -> bool { return freq == FreqSource::HWE; }

    constexpr auto to_byte() const -> uint8_t
    {
        return static_cast<uint8_t>(
            static_cast<unsigned>(scale) | (static_cast<unsigned>(freq) << 1U)
            | (static_cast<unsigned>(proj) << 2U));
    }

    static constexpr auto from_byte(uint8_t b) -> GenotypeProcessMethod
    {
        auto u = static_cast<unsigned>(b);
        return {
            .scale = static_cast<ScaleMethod>(u & 1U),
            .freq = static_cast<FreqSource>((u >> 1U) & 1U),
            .proj = static_cast<Projection>((u >> 2U) & 1U)};
    }

    static constexpr auto StandardizeHWE() -> GenotypeProcessMethod
    {
        return {ScaleMethod::Standardize, FreqSource::HWE, Projection::Direct};
    }
    static constexpr auto CenterHWE() -> GenotypeProcessMethod
    {
        return {ScaleMethod::Center, FreqSource::HWE, Projection::Direct};
    }
    static constexpr auto OrthStandardizeHWE() -> GenotypeProcessMethod
    {
        return {
            ScaleMethod::Standardize, FreqSource::HWE, Projection::Orthogonal};
    }
    static constexpr auto OrthCenterHWE() -> GenotypeProcessMethod
    {
        return {ScaleMethod::Center, FreqSource::HWE, Projection::Orthogonal};
    }
    static constexpr auto Standardize() -> GenotypeProcessMethod
    {
        return {
            ScaleMethod::Standardize, FreqSource::Sample, Projection::Direct};
    }
    static constexpr auto Center() -> GenotypeProcessMethod
    {
        return {ScaleMethod::Center, FreqSource::Sample, Projection::Direct};
    }
    static constexpr auto OrthStandardize() -> GenotypeProcessMethod
    {
        return {
            ScaleMethod::Standardize,
            FreqSource::Sample,
            Projection::Orthogonal};
    }
    static constexpr auto OrthCenter() -> GenotypeProcessMethod
    {
        return {
            ScaleMethod::Center, FreqSource::Sample, Projection::Orthogonal};
    }
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
    static auto format(gelex::GenotypeProcessMethod t, format_context& ctx)
        -> format_context::iterator
    {
        string_view orth = t.is_orthogonal() ? "Orth" : "";
        string_view scale = t.is_center() ? "Center" : "Standardize";
        string_view hwe = t.is_hwe() ? "HWE" : "";
        return fmt::format_to(ctx.out(), "{}{}{}", orth, scale, hwe);
    }
};
}  // namespace fmt

#endif  // GELEX_TYPES_GENOTYPE_PROCESS_METHOD_H_
