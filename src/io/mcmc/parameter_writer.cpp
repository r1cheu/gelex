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

#include "gelex/io/mcmc/parameter_writer.h"

#include <filesystem>

#include "gelex/exception.h"

namespace gelex
{

ParameterWriter::ParameterWriter(
    const mcmc::Result& /*result*/,
    const std::filesystem::path& /*output_path*/)
{
}

auto ParameterWriter::write() -> void
{
    throw GelexException("MCMC parameter writer is legacy and will be removed");
}

}  // namespace gelex
