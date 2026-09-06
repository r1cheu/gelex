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

#ifndef GELEX_BAYES_VARIANCE_DRAWS_H_
#define GELEX_BAYES_VARIANCE_DRAWS_H_

#include <cstdint>
#include <fmt/format.h>
#include <string_view>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/bayes/variance/summary.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

namespace gelex
{

// The total payloads are written even for a single-mode model, where they
// repeat that mode's values, so the layout does not depend on the mode set.
template <GeneticModeSet Modes>
class VarianceSummaryDraws
{
   public:
    VarianceSummaryDraws(BinaryWriter& writer, std::uint64_t draw_count)
        : explained_variance_{
              reserve_per_mode(writer, draw_count, "explained_variance")},
          heritability_{reserve_per_mode(writer, draw_count, "heritability")},
          total_explained_variance_{
              reserve(writer, draw_count, "genetic/total/explained_variance")},
          total_heritability_{
              reserve(writer, draw_count, "genetic/total/heritability")}
    {
    }

    auto append(const VarianceSummary<Modes>& summary) -> void
    {
        explained_variance_.for_each(
            [&]<GeneticMode Mode>(ScalarDraw& draw)
            { draw.append(summary.template genetic<Mode>()); });
        heritability_.for_each(
            [&]<GeneticMode Mode>(ScalarDraw& draw)
            { draw.append(summary.template heritability<Mode>()); });
        total_explained_variance_.append(summary.genetic_total());
        total_heritability_.append(summary.total_heritability());
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto explained_variance() const noexcept -> const ScalarDraw&
    {
        return explained_variance_.template get<Mode>();
    }

    template <GeneticMode Mode>
    [[nodiscard]] auto heritability() const noexcept -> const ScalarDraw&
    {
        return heritability_.template get<Mode>();
    }

    [[nodiscard]] auto total_explained_variance() const noexcept
        -> const ScalarDraw&
    {
        return total_explained_variance_;
    }

    [[nodiscard]] auto total_heritability() const noexcept -> const ScalarDraw&
    {
        return total_heritability_;
    }

   private:
    [[nodiscard]] static auto reserve(
        BinaryWriter& writer,
        std::uint64_t draw_count,
        std::string_view name) -> ScalarDraw
    {
        return ScalarDraw{
            writer.reserve<double>(name, BinaryShape{1, draw_count})};
    }

    [[nodiscard]] static auto reserve_per_mode(
        BinaryWriter& writer,
        std::uint64_t draw_count,
        std::string_view leaf) -> HomogeneousModeValues<Modes, ScalarDraw>
    {
        return generate_mode_values<Modes>(
            [&]<GeneticMode Mode>()
            {
                return reserve(
                    writer,
                    draw_count,
                    fmt::format("genetic/{}/{}", Mode, leaf));
            });
    }

    HomogeneousModeValues<Modes, ScalarDraw> explained_variance_;
    HomogeneousModeValues<Modes, ScalarDraw> heritability_;
    ScalarDraw total_explained_variance_;
    ScalarDraw total_heritability_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_DRAWS_H_
