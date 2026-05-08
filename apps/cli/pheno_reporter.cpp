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

#include "pheno_reporter.h"

#include <fmt/format.h>
#include <string>

#include "cli/report_printer.h"
#include "gelex/infra/logging/formatter.h"

namespace gelex::cli
{

auto PhenoReporter::on_event(const PhenotypeLoadedEvent& event) const -> void
{
    cli::printer().line(
        "   Phenotypes : {} samples ('{}')",
        event.pheno_samples,
        event.trait_name);
}

auto PhenoReporter::on_event(const CovariatesLoadedEvent& event) const -> void
{
    std::string parts;
    if (event.num_quantitative_covariates)
    {
        parts += fmt::format(
            "{} quantitative ({})",
            *event.num_quantitative_covariates,
            gelex::format_names(event.quantitative_names));
    }
    if (event.num_quantitative_covariates && event.num_discrete_covariates)
    {
        parts += ", ";
    }
    if (event.num_discrete_covariates)
    {
        parts += fmt::format(
            "{} discrete ({})",
            *event.num_discrete_covariates,
            gelex::format_names(event.discrete_names));
    }
    cli::printer().line("   Covariates : {}", parts);
}

}  // namespace gelex::cli
