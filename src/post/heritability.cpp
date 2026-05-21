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

#include "gelex/post/heritability.h"

#include <string>

#include <fmt/format.h>
#include <Eigen/Core>

#include "gelex/infra/stats/diagnostics.h"
#include "gelex/io/detail/binary_reader.h"
#include "gelex/model/bayes/legacy_algorithm_shape.h"
#include "gelex/post/detail/utils.h"
#include "gelex/types/genetic_effect_type.h"

namespace gelex
{

HeritabilityPosteriorProcessor::HeritabilityPosteriorProcessor(
    std::span<const io::detail::BinaryReader> readers,
    double hdpi_threshold,
    const stats::Chains& genetic_variances,
    std::span<const GeneticMode> kinds)
    : readers_{readers},
      hdpi_threshold_{hdpi_threshold},
      genetic_variances_{genetic_variances},
      kinds_{kinds}
{
}

auto HeritabilityPosteriorProcessor::process() -> std::vector<ParameterDiag>
{
    if (genetic_variances_.empty())
    {
        return {};
    }

    const auto n_chains = readers_.size();
    const auto n_kinds = genetic_variances_.front().rows() - 1;

    auto phenotypic_var = assemble_phenotypic_variance();

    const auto n_rows = genetic_variances_.front().rows();
    stats::Chains h2_chains(
        n_chains, Eigen::MatrixXd(n_rows, phenotypic_var[0].cols()));
    for (size_t ci = 0; ci < n_chains; ++ci)
    {
        h2_chains[ci] = (genetic_variances_[ci].array().rowwise()
                         / phenotypic_var[ci].array().row(0))
                            .matrix();
    }

    auto diags
        = post::detail::compute_posterior_summaries(h2_chains, hdpi_threshold_);

    for (Eigen::Index ki = 0; ki < n_kinds; ++ki)
    {
        auto& diag = diags[static_cast<size_t>(ki)];
        diag.section = "Genetic";
        diag.name = std::string(
            bayes::to_heritability_label(kinds_[static_cast<size_t>(ki)]));
    }
    diags.back().section = "Genetic";
    diags.back().name = "H²";

    return diags;
}

auto HeritabilityPosteriorProcessor::assemble_phenotypic_variance() const
    -> stats::Chains
{
    const auto n_chains = readers_.size();
    const auto total_row = genetic_variances_.front().rows() - 1;

    stats::Chains pheno_var;
    pheno_var.reserve(n_chains);
    for (size_t ci = 0; ci < n_chains; ++ci)
    {
        pheno_var.emplace_back(genetic_variances_[ci].row(total_row));
    }

    auto resid_path = fmt::format("{}/variance", EffectType::residual());
    for (size_t ci = 0; ci < n_chains; ++ci)
    {
        pheno_var[ci] += readers_[ci].to_map<double>(resid_path);
    }

    auto paths = readers_.front().section_paths();
    auto random_var_prefix = fmt::format("{}/variance/", EffectType::random());
    for (auto path : paths)
    {
        if (!path.starts_with(random_var_prefix))
        {
            continue;
        }
        for (size_t ci = 0; ci < n_chains; ++ci)
        {
            pheno_var[ci] += readers_[ci].to_map<double>(path);
        }
    }

    return pheno_var;
}

}  // namespace gelex
