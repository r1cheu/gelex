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

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <CLI/CLI.hpp>
#include <Eigen/Core>

#include "gelex/data/covariates.h"
#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/index.h"
#include "gelex/data/reader.h"
#include "gelex/types/fixed_designs.h"

namespace cli
{
struct BaseDataConfig
{
    std::string pheno_path;
    int pheno_col{0};
    std::optional<std::string> qcovar_path;
    std::optional<std::string> dcovar_path;
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
    std::vector<gelex::dataframe::Index<std::string>*>& indices,
    const gelex::dataframe::Index<std::string>& common_index) {
    handler.load_indices(indices);
    handler.gather(common_index);
    std::move(handler).results();
};

template <BaseDataHandler Handler>
auto load_base_data(Handler& handler, const BaseDataConfig& config) -> BaseData
{
    std::vector<gelex::dataframe::Index<std::string>*> indices;

    auto pheno_col_offset = static_cast<std::size_t>(config.pheno_col);
    auto phenotype = gelex::read_pheno(config.pheno_path, &pheno_col_offset);
    indices.push_back(&phenotype.index());

    std::optional<gelex::dataframe::DataFrame<std::string>> qcovar;
    std::optional<gelex::dataframe::DataFrame<std::string>> dcovar;
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

    auto common_index = gelex::dataframe::intersect<std::string>(indices);

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
    return BaseData{
        .phenotype = phenotype.col(0).to_mat<double>(),
        .fixed_design = std::move(*fixed_design),
        .sample_ids = std::move(common_index).take_keys(),
    };
}

template <BaseDataHandler Handler>
auto load_base_data(Handler& handler, CLI::App& cmd) -> BaseData
{
    BaseDataConfig config{
        .pheno_path = cmd.get_option("--pheno")->as<std::string>(),
        .pheno_col = cmd.get_option("--pheno-col")->as<int>(),
        .qcovar_path = cmd.get_option("--qcovar")->count() > 0
                           ? std::make_optional(
                                 cmd.get_option("--qcovar")->as<std::string>())
                           : std::nullopt,
        .dcovar_path = cmd.get_option("--dcovar")->count() > 0
                           ? std::make_optional(
                                 cmd.get_option("--dcovar")->as<std::string>())
                           : std::nullopt,
    };
    return load_base_data(handler, config);
}

}  // namespace cli

#endif  // APPS_CLI_COMMON_DATA_H_
