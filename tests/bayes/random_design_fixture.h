// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_TEST_BAYES_RANDOM_DESIGN_FIXTURE_H_
#define GELEX_TEST_BAYES_RANDOM_DESIGN_FIXTURE_H_

#include <Eigen/Core>
#include <fmt/format.h>
#include <iterator>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "gelex/bayes/design.h"
#include "gelex/data/dataframe/reader.h"

#include "file_fixture.h"

namespace gelex::test
{

inline auto make_random_design(
    std::string name,
    std::span<const std::string> column_names,
    const Eigen::Ref<const Eigen::MatrixXd>& X) -> bayes::RandomDesign
{
    std::string content = "id";
    for (const auto& column_name : column_names)
    {
        fmt::format_to(std::back_inserter(content), "\t{}", column_name);
    }
    content.push_back('\n');

    for (Eigen::Index row = 0; row < X.rows(); ++row)
    {
        fmt::format_to(std::back_inserter(content), "s{}", row);
        for (Eigen::Index column = 0; column < X.cols(); ++column)
        {
            fmt::format_to(std::back_inserter(content), "\t{}", X(row, column));
        }
        content.push_back('\n');
    }

    FileFixture files;
    const auto path = files.create_text_file(content, ".tsv");
    ReadOptions options;
    options.index_cols = {0};
    const auto frame
        = read_dataframe<std::string, double>(path.string(), options);
    return bayes::make_quantitative_random_design(frame, std::move(name));
}

}  // namespace gelex::test

#endif  // GELEX_TEST_BAYES_RANDOM_DESIGN_FIXTURE_H_
