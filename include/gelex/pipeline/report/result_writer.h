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

#ifndef GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_
#define GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_

#include <filesystem>

namespace gelex
{

class MCMCResult;

class MCMCResultWriter
{
   public:
    MCMCResultWriter(
        const MCMCResult& result,
        const std::filesystem::path& bim_file_path);

    auto save(const std::filesystem::path& prefix) const -> void;

   private:
    const MCMCResult* result_;
    std::filesystem::path bim_file_path_;
};

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_BAYES_RESULT_WRITER_H_
