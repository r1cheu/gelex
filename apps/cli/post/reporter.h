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

#ifndef APPS_CLI_POST_REPORTER_H_
#define APPS_CLI_POST_REPORTER_H_

#include <string>
#include <vector>

#include "gelex/post/diagnostic.h"

namespace cli
{

class PostReporter
{
   public:
    auto show_banner() const -> void;
    auto show_start(const std::vector<std::string>& in_prefixes) const -> void;
    auto show_diagnostics(
        const std::vector<gelex::ParameterDiag>& diags,
        double hdpi_prob) const -> void;
};

}  // namespace cli

#endif  // APPS_CLI_POST_REPORTER_H_
