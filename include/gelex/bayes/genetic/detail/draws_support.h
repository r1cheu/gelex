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

#ifndef GELEX_BAYES_GENETIC_DETAIL_DRAWS_SUPPORT_H_
#define GELEX_BAYES_GENETIC_DETAIL_DRAWS_SUPPORT_H_

#include <Eigen/Core>
#include <cstddef>
#include <cstdint>
#include <fmt/format.h>
#include <string>
#include <string_view>
#include <utility>

#include "gelex/bayes/basic_draw.h"
#include "gelex/bayes/genetic/traits.h"
#include "gelex/bayes/genetic_family.h"
#include "gelex/io/binary_format.h"
#include "gelex/io/binary_writer.h"

namespace gelex::detail
{

struct GeneticDrawsDimensions
{
    Eigen::Index marker_count;
    std::uint64_t draw_count;
};

class GeneticDrawsBuilder
{
   public:
    GeneticDrawsBuilder(
        BinaryWriter& writer,
        std::string prefix,
        GeneticDrawsDimensions dimensions)
        : writer_{&writer}, prefix_{std::move(prefix)}, dimensions_{dimensions}
    {
    }

    [[nodiscard]] auto scalar(std::string_view name) -> ScalarDraw
    {
        return ScalarDraw{reserve<double>(name, 1)};
    }

    [[nodiscard]] auto vector(std::string_view name, Eigen::Index rows)
        -> VectorDraw
    {
        return VectorDraw{reserve<float>(name, rows)};
    }

    template <std::size_t CategoryCount>
    [[nodiscard]] auto category(std::string_view name, Eigen::Index rows)
        -> CategoryDraw<CategoryCount>
    {
        return CategoryDraw<CategoryCount>{reserve<std::uint8_t>(name, rows)};
    }

    [[nodiscard]] auto marker_count() const noexcept -> Eigen::Index
    {
        return dimensions_.marker_count;
    }

   private:
    template <SupportedDtype T>
    [[nodiscard]] auto reserve(std::string_view name, Eigen::Index rows)
        -> PayloadWriter<T>
    {
        return writer_->reserve<T>(
            fmt::format("{}/{}", prefix_, name),
            BinaryShape{
                static_cast<std::uint64_t>(rows), dimensions_.draw_count});
    }

    BinaryWriter* writer_;
    std::string prefix_;
    GeneticDrawsDimensions dimensions_;
};

template <VarianceLayout Kind>
[[nodiscard]] auto make_marker_variance_draw(GeneticDrawsBuilder& builder)
    -> marker_variance_draw_t<Kind>
{
    if constexpr (Kind == VarianceLayout::Pooled)
    {
        return builder.scalar("variance");
    }
    else
    {
        return builder.vector("variance", builder.marker_count());
    }
}

template <MixtureWeightUpdate Update>
[[nodiscard]] auto make_probability_draw(
    GeneticDrawsBuilder& builder,
    std::string_view name) -> weight_draw_t<Update, ScalarDraw>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
    {
        return builder.scalar(name);
    }
    else
    {
        return EmptyDraw{};
    }
}

template <MixtureWeightUpdate Update, std::size_t ClassCount>
[[nodiscard]] auto make_probabilities_draw(GeneticDrawsBuilder& builder)
    -> weight_draw_t<Update, VectorDraw>
{
    if constexpr (Update == MixtureWeightUpdate::Enabled)
    {
        return builder.vector(
            "probabilities", static_cast<Eigen::Index>(ClassCount));
    }
    else
    {
        return EmptyDraw{};
    }
}
// rows non-additive: they are shares, not a decomposition of genetic variance.
template <std::size_t ComponentCount>
[[nodiscard]] auto make_component_explained_variance_draw(
    GeneticDrawsBuilder& builder) -> VectorDraw
{
    return builder.vector(
        "component_explained_variance",
        static_cast<Eigen::Index>(ComponentCount));
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_GENETIC_DETAIL_DRAWS_SUPPORT_H_
