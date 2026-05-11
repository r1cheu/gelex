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

#include "gelex/model/bayes/algorithm_shape.h"

#include <algorithm>
#include <optional>
#include <span>
#include <string_view>

#include "gelex/exception.h"
#include "gelex/model/bayes/bayes_policy.h"

namespace gelex::bayes
{

auto try_resolve_shape(
    const BayesPolicy& policy,
    std::span<const GeneticMode> requested) -> std::optional<AlgorithmShape>
{
    const bool has_a = std::ranges::contains(requested, GeneticMode::A);
    const bool has_d = std::ranges::contains(requested, GeneticMode::D);

    auto matches_request = [&](AlgorithmShape s)
    {
        if (has_a && has_d)
        {
            return s == AlgorithmShape::ad_joint
                   || s == AlgorithmShape::ad_independent;
        }
        return s == (has_a ? AlgorithmShape::a_only : AlgorithmShape::d_only);
    };

    const auto it = std::ranges::find_if(policy.shapes, matches_request);
    if (it == policy.shapes.end())
    {
        return std::nullopt;
    }
    return *it;
}

auto resolve_shape(
    const BayesPolicy& policy,
    std::span<const GeneticMode> requested) -> AlgorithmShape
{
    if (auto shape = try_resolve_shape(policy, requested))
    {
        return *shape;
    }
    throw GelexException{
        "resolve_shape: policy does not support the requested effect set"};
}

auto to_variance_label(AlgorithmShape shape) -> std::string_view
{
    switch (shape)
    {
        case AlgorithmShape::a_only:
            return "σ²_add";
        case AlgorithmShape::d_only:
            return "σ²_dom";
        case AlgorithmShape::ad_independent:
        case AlgorithmShape::ad_joint:
            return "σ²_g";
    }
    return "unknown";
}

auto to_heritability_label(AlgorithmShape shape) -> std::string_view
{
    switch (shape)
    {
        case AlgorithmShape::a_only:
            return "h²";
        case AlgorithmShape::d_only:
            return "δ²";
        case AlgorithmShape::ad_independent:
        case AlgorithmShape::ad_joint:
            return "H²";
    }
    return "unknown";
}

auto to_file_suffix(AlgorithmShape shape) -> std::string_view
{
    switch (shape)
    {
        case AlgorithmShape::a_only:
            return "add";
        case AlgorithmShape::d_only:
            return "dom";
        case AlgorithmShape::ad_independent:
        case AlgorithmShape::ad_joint:
            return "ad";
    }
    return "unknown";
}

auto to_variance_label(GeneticMode mode) -> std::string_view
{
    switch (mode)
    {
        case GeneticMode::A:
            return "σ²_add";
        case GeneticMode::D:
            return "σ²_dom";
    }
    return "unknown";
}

auto to_heritability_label(GeneticMode mode) -> std::string_view
{
    switch (mode)
    {
        case GeneticMode::A:
            return "h²";
        case GeneticMode::D:
            return "δ²";
    }
    return "unknown";
}

auto to_file_suffix(GeneticMode mode) -> std::string_view
{
    switch (mode)
    {
        case GeneticMode::A:
            return "add";
        case GeneticMode::D:
            return "dom";
    }
    return "unknown";
}

}  // namespace gelex::bayes
