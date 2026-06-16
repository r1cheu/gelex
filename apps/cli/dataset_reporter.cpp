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

#include "dataset_reporter.h"

#include <cstddef>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"

namespace cli
{

auto DatasetReporter::show_section() const -> void
{
    cli::printer().block(gelex::section("[Dataset Summary]"));
}

auto DatasetReporter::show_intersection(size_t common_samples) const -> void
{
    cli::printer().line("   Intersection : {} common samples", common_samples);
}

}  // namespace cli
