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

#include "gelex/model/bayes/writer/covariate_writer.h"

#include <fmt/format.h>
#include <memory>

#include "gelex/io/detail/text_writer.h"
#include "gelex/types/mcmc_result.h"

namespace gelex
{

CovariateWriter::CovariateWriter(
    const FixedSummary& fixed,
    std::span<const RandomSummary> random,
    const std::filesystem::path& output_path)
    : fixed_(&fixed),
      random_(random),
      writer_(std::make_unique<io::detail::TextWriter>(output_path))
{
}

CovariateWriter::~CovariateWriter() = default;

auto CovariateWriter::write() -> void
{
    writer_->write_header({"term", "mean", "stddev"});

    write_fixed_effects();
    write_random_effects();
}

auto CovariateWriter::write_fixed_effects() -> void
{
    fixed_->for_each_term(
        [&](const std::string& term, Eigen::Index i)
        {
            writer_->write(
                fmt::format(
                    "{}\t{}\t{}",
                    term,
                    fixed_->coeffs.mean(i),
                    fixed_->coeffs.stddev(i)));
        });
}

auto CovariateWriter::write_random_effects() -> void
{
    for (const auto& rand : random_)
    {
        rand.for_each_term(
            [&](const std::string& term, Eigen::Index i)
            {
                writer_->write(
                    fmt::format(
                        "{}\t{}\t{}",
                        term,
                        rand.coeffs.mean(i),
                        rand.coeffs.stddev(i)));
            });
    }
}

}  // namespace gelex
