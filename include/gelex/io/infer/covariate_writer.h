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

#ifndef GELEX_IO_INFER_COVARIATE_WRITER_H_
#define GELEX_IO_INFER_COVARIATE_WRITER_H_

#include <filesystem>
#include <memory>
#include <span>

namespace gelex::io::detail
{
class TextWriter;
}

namespace gelex
{

struct FixedSummary;
struct RandomSummary;

class CovariateWriter
{
   public:
    CovariateWriter(
        const FixedSummary& fixed,
        std::span<const RandomSummary> random,
        const std::filesystem::path& output_path);
    ~CovariateWriter();
    CovariateWriter(const CovariateWriter&) = delete;
    CovariateWriter& operator=(const CovariateWriter&) = delete;
    CovariateWriter(CovariateWriter&&) = delete;
    CovariateWriter& operator=(CovariateWriter&&) = delete;

    auto write() -> void;

   private:
    auto write_fixed_effects() -> void;
    auto write_random_effects() -> void;

    const FixedSummary* fixed_;
    std::span<const RandomSummary> random_;
    std::unique_ptr<io::detail::TextWriter> writer_;
};

}  // namespace gelex

#endif  // GELEX_IO_INFER_COVARIATE_WRITER_H_
