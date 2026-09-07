// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_BAYES_VARIANCE_DRAWS_H_
#define GELEX_BAYES_VARIANCE_DRAWS_H_

#include <cstdint>
#include <fmt/format.h>
#include <string_view>

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
    using payload_type = PayloadWriter<double>;

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
            [&]<GeneticMode Mode>(payload_type& draw)
            { draw.append(summary.template genetic<Mode>()); });
        heritability_.for_each(
            [&]<GeneticMode Mode>(payload_type& draw)
            { draw.append(summary.template heritability<Mode>()); });
        total_explained_variance_.append(summary.genetic_total());
        total_heritability_.append(summary.total_heritability());
    }

   private:
    [[nodiscard]] static auto reserve(
        BinaryWriter& writer,
        std::uint64_t draw_count,
        std::string_view name) -> payload_type
    {
        return writer.reserve<double>(name, BinaryShape{1, draw_count});
    }

    [[nodiscard]] static auto reserve_per_mode(
        BinaryWriter& writer,
        std::uint64_t draw_count,
        std::string_view leaf) -> HomogeneousModeValues<Modes, payload_type>
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

    HomogeneousModeValues<Modes, payload_type> explained_variance_;
    HomogeneousModeValues<Modes, payload_type> heritability_;
    payload_type total_explained_variance_;
    payload_type total_heritability_;
};

}  // namespace gelex

#endif  // GELEX_BAYES_VARIANCE_DRAWS_H_
