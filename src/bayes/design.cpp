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

#include "gelex/bayes/design.h"

#include <array>
#include <span>
#include <string>
#include <vector>

#include <fmt/format.h>

#include "gelex/infra/field_flag.h"
#include "gelex/infra/field_visitor.h"

namespace gelex::bayes
{

auto RandomDesign::visit(infra::FieldVisitor& visitor) const -> void
{
    std::vector<std::string> coeff_names;
    coeff_names.reserve(levels.size());
    for (const auto& level : levels)
    {
        coeff_names.push_back(name.empty() ? level : name + "_" + level);
    }

    visitor.on(
        "coeffs_names",
        std::span<const std::string>{coeff_names},
        FieldFlag::summary);

    const std::array<std::string, 1> variance_names{fmt::format("σ²_{}", name)};
    visitor.on(
        "variance_names",
        std::span<const std::string>{variance_names},
        FieldFlag::summary);
}

}  // namespace gelex::bayes
