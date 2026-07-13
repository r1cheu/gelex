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

#ifndef APPS_CLI_COMMON_DATA_H_
#define APPS_CLI_COMMON_DATA_H_

#include <Eigen/Core>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/exception.h"
#include "gelex/infra/stats/rank_inverse_norm_transform.h"
#include "gelex/types/fixed_designs.h"

#include "cli/logging.h"

namespace cli
{
struct BaseDataConfig
{
    std::string pheno_path;
    int pheno_col{0};
    std::optional<std::string> qcovar_path;
    std::optional<std::string> dcovar_path;
    gelex::RintType transform{gelex::RintType::None};
    double int_offset{3.0 / 8.0};
};

struct BaseData
{
    Eigen::VectorXd phenotype;
    gelex::FixedDesign fixed_design;
    std::vector<std::string> sample_ids;
};

template <typename T>
concept BaseDataHandler = requires(
    T& handler,
    std::vector<const gelex::DataFrameIndex<std::string>*>& indices,
    const gelex::DataFrameIndex<std::string>& common_index) {
    handler.load_indices(indices);
    handler.gather(common_index);
    std::move(handler).results();
};

template <BaseDataHandler Handler>
auto load_base_data(Handler& handler, const BaseDataConfig& config) -> BaseData
{
    std::vector<const gelex::DataFrameIndex<std::string>*> indices;

    auto pheno_col_offset = static_cast<std::size_t>(config.pheno_col);
    auto phenotype = gelex::read_pheno(config.pheno_path, &pheno_col_offset);
    indices.push_back(&phenotype.index());

    std::optional<gelex::DataFrame<std::string>> qcovar;
    std::optional<gelex::DataFrame<std::string>> dcovar;
    if (config.qcovar_path)
    {
        qcovar = std::make_optional(gelex::read_qcovar(*config.qcovar_path));
        indices.push_back(&qcovar->index());
    }
    if (config.dcovar_path)
    {
        dcovar = std::make_optional(gelex::read_dcovar(*config.dcovar_path));
        indices.push_back(&dcovar->index());
    }
    handler.load_indices(indices);

    auto common_index = gelex::intersect<std::string>(indices);
    if (common_index.size() == 0)
    {
        throw gelex::GelexException(
            "No samples remain after intersecting phenotype, covariate, and "
            "GRM IDs; check that FID/IID match across all input files");
    }

    phenotype.gather(common_index);
    if (qcovar)
    {
        qcovar->gather(common_index);
    }
    if (dcovar)
    {
        dcovar->gather(common_index);
    }
    handler.gather(common_index);

    std::optional<gelex::FixedDesign> fixed_design;

    if (!qcovar && !dcovar)
    {
        fixed_design = std::make_optional(
            gelex::FixedDesign::make(
                static_cast<Eigen::Index>(common_index.size())));
    }
    else
    {
        fixed_design = std::make_optional(
            gelex::FixedDesign::make(
                qcovar ? std::make_optional(
                             gelex::make_quantitative_covariate(*qcovar))
                       : std::nullopt,
                dcovar ? std::make_optional(
                             gelex::make_discrete_covariate(*dcovar))
                       : std::nullopt));
    }
    auto pheno_vec = phenotype.col(0).to_mat<double>();

    switch (config.transform)
    {
        case gelex::RintType::None:
            break;
        case gelex::RintType::Direct:
            cli::logging::get()->info(
                "   Method: Direct INT (DINT), offset (k): {}",
                config.int_offset);
            gelex::direct_int(pheno_vec, config.int_offset);
            break;
        case gelex::RintType::Indirect:
            cli::logging::get()->info(
                "   Method: Indirect INT (IINT), offset (k): {}",
                config.int_offset);
            gelex::indirect_int(pheno_vec, fixed_design->X, config.int_offset);
            fixed_design = gelex::FixedDesign::make(pheno_vec.size());
            break;
    }

    return BaseData{
        .phenotype = std::move(pheno_vec),
        .fixed_design = std::move(*fixed_design),
        .sample_ids = std::move(common_index).take_keys(),
    };
}

}  // namespace cli

#endif  // APPS_CLI_COMMON_DATA_H_
