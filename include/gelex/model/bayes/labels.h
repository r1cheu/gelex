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

#ifndef GELEX_MODEL_BAYES_LABELS_H_
#define GELEX_MODEL_BAYES_LABELS_H_

#include <string_view>

#include "gelex/types/genetic_effect_type.h"

namespace gelex::bayes
{

inline auto to_variance_label(GeneticMode mode) -> std::string_view
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

inline auto to_heritability_label(GeneticMode mode) -> std::string_view
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

inline auto to_file_suffix(GeneticMode mode) -> std::string_view
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

#endif  // GELEX_MODEL_BAYES_LABELS_H_
