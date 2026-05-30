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

#include "gelex/bayes/recipe_options.h"

#include <optional>
#include <utility>

#include <fmt/format.h>

#include "gelex/exception.h"
#include "gelex/types/constrained_value.h"
#include "gelex/types/constrained_vector.h"

namespace gelex::bayes
{

EffectConfig::EffectConfig(
    std::optional<OpenUnitInterval<double>> heritability,
    std::optional<Simplex<double>> proportion,
    std::optional<ScaleMultiplier<double>> multiplier,
    std::optional<bool> proportion_update)
    : heritability_(heritability),
      proportion_(std::move(proportion)),
      multiplier_(std::move(multiplier)),
      proportion_update_(proportion_update)
{
    if (proportion_ && multiplier_
        && proportion_->size() != multiplier_->size())
    {
        throw GelexException(
            fmt::format(
                "EffectConfig: proportion size {} != multiplier size {}",
                proportion_->size(),
                multiplier_->size()));
    }
}

}  // namespace gelex::bayes
