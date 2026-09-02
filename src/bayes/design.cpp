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

#include "gelex/bayes/design.h"

#include <Eigen/Core>
#include <cstddef>
#include <fmt/format.h>
#include <string>
#include <utility>
#include <vector>

#include "gelex/data/dataframe/dataframe.h"
#include "gelex/data/dataframe/encode.h"
#include "gelex/exception.h"

namespace gelex::bayes
{

RandomDesign::RandomDesign(
    std::string name,
    std::vector<std::string> column_names,
    Eigen::MatrixXd X)
    : name_(std::move(name)),
      column_names_(std::move(column_names)),
      matrix_(std::move(X)),
      xtx_diag_(matrix_.colwise().squaredNorm().transpose())
{
    if (name_.empty())
    {
        throw GelexException("RandomDesign: name must not be empty");
    }
    if (matrix_.rows() == 0)
    {
        throw GelexException(
            fmt::format("RandomDesign '{}': X must not have zero rows", name_));
    }
    if (matrix_.cols() == 0)
    {
        throw GelexException(
            fmt::format(
                "RandomDesign '{}': X must not have zero columns", name_));
    }
    if (!xtx_diag_.allFinite() || (xtx_diag_.array() <= 0.0).any())
    {
        throw GelexException(
            fmt::format(
                "RandomDesign '{}': X column squared norms must be finite and "
                "positive",
                name_));
    }
}

auto make_random_designs(const DataFrame<std::string>& frame)
    -> std::vector<RandomDesign>
{
    std::vector<RandomDesign> random_designs;
    random_designs.reserve(frame.cols());
    for (std::size_t i = 0; i < frame.cols(); ++i)
    {
        auto result = one_hot_encode(frame.col(i));
        random_designs.push_back(
            RandomDesign{
                std::move(result.name),
                std::move(result.level_names),
                std::move(result.data)});
    }
    return random_designs;
}

auto make_quantitative_random_design(
    const DataFrame<std::string>& frame,
    std::string name) -> RandomDesign
{
    return RandomDesign{std::move(name), frame.names(), frame.to_mat<double>()};
}

}  // namespace gelex::bayes
