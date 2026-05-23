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

#include "gelex/post/residual.h"

#include <string_view>

#include "gelex/io/detail/binary_reader.h"
#include "gelex/post/detail/utils.h"

namespace gelex
{

ResidualPosteriorProcessor::ResidualPosteriorProcessor(
    std::span<const io::detail::BinaryReader> readers,
    double hdpi_threshold)
    : readers_{readers}, hdpi_threshold_{hdpi_threshold}
{
}

auto ResidualPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    constexpr std::string_view resid_path = "residual/0/variance";
    return {post::detail::summarize_section(
        readers_, resid_path, hdpi_threshold_, "Residual", "σ²_e")};
}

}  // namespace gelex
