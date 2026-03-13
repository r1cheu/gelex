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

#include "gelex/predict/covariate_reader.h"

#include <format>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gelex/exception.h"
#include "gelex/io/parser.h"

namespace gelex
{

auto CovariateReader::parse(const std::filesystem::path& path)
    -> CovariateParams
{
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);

    std::string line;
    if (!std::getline(file, line))
    {
        throw FileFormatException(
            std::format(
                "parameter file '{}' is empty or has no header",
                path.string()));
    }

    std::vector<std::string> term_names;
    std::vector<double> coefficients;

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::vector<std::string> fields;
        std::istringstream row_stream(line);
        for (std::string token; std::getline(row_stream, token, '\t');)
        {
            fields.emplace_back(token);
        }

        if (fields.size() < 2 || fields[0].empty() || fields[1].empty())
        {
            continue;
        }

        double mean{};
        try
        {
            mean = std::stod(fields[1]);
        }
        catch (const std::exception&)
        {
            continue;
        }

        term_names.push_back(fields[0]);
        coefficients.push_back(mean);
    }

    if (term_names.empty())
    {
        throw FileFormatException(
            std::format(
                "no valid data rows in parameter file '{}'", path.string()));
    }

    CovariateParams result;
    result.term_names = std::move(term_names);
    result.coefficients = Eigen::Map<Eigen::VectorXd>(
        coefficients.data(), static_cast<Eigen::Index>(coefficients.size()));

    return result;
}

}  // namespace gelex
