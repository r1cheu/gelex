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

#ifndef GELEX_TEST_PREDICT_ENGINE_FIXTURE_H
#define GELEX_TEST_PREDICT_ENGINE_FIXTURE_H

#include <fstream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "bed_fixture.h"
#include "file_fixture.h"

namespace gelex::test
{

class PredictEngineTestFixture : public BedFixture
{
   public:
    explicit PredictEngineTestFixture() = default;

    static std::string write_row(std::span<const std::string> row)
    {
        std::string line;
        for (size_t i = 0; i < row.size(); ++i)
        {
            line += row[i];
            if (i < row.size() - 1)
            {
                line += "\t";
            }
        }
        line += "\n";
        return line;
    }

    static std::string create_table(
        const std::vector<std::string>& headers,
        const std::vector<std::vector<std::string>>& rows)
    {
        std::string content = write_row(headers);
        for (const auto& row : rows)
        {
            content += write_row(row);
        }
        return content;
    }

    std::filesystem::path create_snp_effects_file(
        const std::vector<std::vector<std::string>>& snp_rows,
        bool has_dominance = false)
    {
        std::vector<std::string> headers
            = {"Chrom", "Position", "ID", "A1", "A2", "A1Freq", "Add"};
        if (has_dominance)
        {
            headers.emplace_back("Dom");
        }
        auto content = create_table(headers, snp_rows);
        return get_file_fixture().create_text_file(content, ".snp.eff");
    }

    std::filesystem::path create_param_file(
        const std::vector<std::vector<std::string>>& rows)
    {
        return get_file_fixture().create_text_file(
            create_table(
                {"term",
                 "mean",
                 "stddev",
                 "percentile_5",
                 "percentile_95",
                 "ess",
                 "rhat"},
                rows),
            ".param");
    }

    std::filesystem::path create_param_intercept_only(double intercept)
    {
        std::vector<std::vector<std::string>> rows
            = {{"Intercept",
                std::to_string(intercept),
                "0.1",
                "0.8",
                "1.2",
                "1000",
                "1.0"}};
        return create_param_file(rows);
    }

    std::filesystem::path create_param_with_qcovar(
        double intercept,
        const std::vector<std::pair<std::string, double>>& coefs)
    {
        std::vector<std::vector<std::string>> rows;
        rows.push_back(
            {"Intercept",
             std::to_string(intercept),
             "0.1",
             "0.8",
             "1.2",
             "1000",
             "1.0"});
        for (const auto& [name, coef] : coefs)
        {
            rows.push_back(
                {name,
                 std::to_string(coef),
                 "0.05",
                 "0.1",
                 "0.3",
                 "800",
                 "1.01"});
        }
        return create_param_file(rows);
    }

    std::filesystem::path create_param_with_dcovar(
        double intercept,
        const std::vector<std::pair<std::string, double>>& coefs)
    {
        std::vector<std::vector<std::string>> rows;
        rows.push_back(
            {"Intercept",
             std::to_string(intercept),
             "0.1",
             "0.8",
             "1.2",
             "1000",
             "1.0"});
        for (const auto& [name, coef] : coefs)
        {
            rows.push_back(
                {name,
                 std::to_string(coef),
                 "0.02",
                 "-0.34",
                 "-0.26",
                 "900",
                 "1.02"});
        }
        return create_param_file(rows);
    }

    std::filesystem::path create_param_full(
        double intercept,
        const std::vector<std::pair<std::string, double>>& qcovar_coefs,
        const std::vector<std::pair<std::string, double>>& dcovar_coefs)
    {
        std::vector<std::vector<std::string>> rows;
        rows.push_back(
            {"Intercept",
             std::to_string(intercept),
             "0.1",
             "0.8",
             "1.2",
             "1000",
             "1.0"});
        for (const auto& [name, coef] : qcovar_coefs)
        {
            rows.push_back(
                {name,
                 std::to_string(coef),
                 "0.05",
                 "0.1",
                 "0.3",
                 "800",
                 "1.01"});
        }
        for (const auto& [name, coef] : dcovar_coefs)
        {
            rows.push_back(
                {name,
                 std::to_string(coef),
                 "0.02",
                 "-0.34",
                 "-0.26",
                 "900",
                 "1.02"});
        }
        return create_param_file(rows);
    }

    template <typename T>
    std::filesystem::path create_covar_file(
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        const std::vector<std::pair<std::string, std::vector<T>>>& covars,
        std::string_view suffix)
    {
        std::string content = "FID\tIID";
        for (const auto& cov : covars)
        {
            content += "\t" + cov.first;
        }
        content += "\n";
        for (size_t i = 0; i < iids.size(); ++i)
        {
            content += fids[i] + "\t" + iids[i];
            for (const auto& cov : covars)
            {
                content += "\t";
                if constexpr (std::is_same_v<T, double>)
                {
                    content += std::format("{}", cov.second[i]);
                }
                else
                {
                    content += cov.second[i];
                }
            }
            content += "\n";
        }
        return get_file_fixture().create_text_file(content, suffix);
    }

    std::filesystem::path create_qcovar_file(
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        const std::vector<std::pair<std::string, std::vector<double>>>& qcovars)
    {
        return create_covar_file(fids, iids, qcovars, ".qcovar");
    }

    std::filesystem::path create_dcovar_file(
        const std::vector<std::string>& fids,
        const std::vector<std::string>& iids,
        const std::vector<std::pair<std::string, std::vector<std::string>>>&
            dcovars)
    {
        return create_covar_file(fids, iids, dcovars, ".dcovar");
    }

    static std::pair<std::vector<std::string>, std::vector<std::string>>
    read_fam(const std::filesystem::path& path)
    {
        std::vector<std::string> fids;
        std::vector<std::string> iids;
        std::ifstream fam_file(path);
        std::string line;
        while (std::getline(fam_file, line))
        {
            std::istringstream iss(line);
            std::string buffer;
            iss >> buffer;
            fids.push_back(buffer);
            iss >> buffer;
            iids.push_back(buffer);
        }
        return {std::move(fids), std::move(iids)};
    }
};

}  // namespace gelex::test

#endif  // GELEX_TEST_PREDICT_ENGINE_FIXTURE_H
