#include "bim_loader.h"

#include <filesystem>
#include <fstream>
#include <ranges>
#include <string_view>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"

#include "../src/data/parser.h"

namespace gelex::detail
{

BimLoader::BimLoader(const std::filesystem::path& path)
{
    auto file = detail::open_file<std::ifstream>(path, std::ios_base::in);
    char delimiter = detect_delimiter(file);
    try
    {
        set_snp_info(delimiter, file);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}

void BimLoader::set_snp_info(char delimiter, std::ifstream& file)
{
    snp_info_.clear();
    snp_info_.reserve(1024);  // Start small
    std::string line;
    int n_line = 0;
    std::array<std::string_view, 6> cols;

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }

        auto tokens = line | std::views::split(delimiter)
                      | std::views::transform([](auto&& rng)
                                              { return std::string_view(rng); })
                      | std::views::filter([](std::string_view sv)
                                           { return !sv.empty(); })
                      | std::views::take(6);

        auto [in, out] = std::ranges::copy(tokens, cols.begin());

        if (out != cols.end())
        {
            const auto count = std::distance(cols.begin(), out);
            throw InconsistentColumnCountException(
                std::format(
                    "Line {}: has {} columns, expected 6", n_line, count));
        }
        try
        {
            snp_info_.emplace_back(
                SnpInfo{
                    .chrom = cols[0][0],
                    .id = std::string(cols[1]),
                    .base_coordinate = detail::parse_number<int>(cols[3]),
                    .A1 = cols[4][0],
                    .A2 = cols[5][0],
                });
        }
        catch (const gelex::GelexException& err)
        {
            throw DataParseException(std::format("{}: {}", n_line, err.what()));
        };
    }
}

char BimLoader::detect_delimiter(std::ifstream& file)
{
    std::string probe_line;

    // escape empty lines
    while (std::getline(file, probe_line) && probe_line.empty())
    {
    }
    bool is_tab = !probe_line.empty() && probe_line.contains('\t');

    file.clear();
    file.seekg(0);

    return is_tab ? '\t' : ' ';
}

std::vector<std::string> BimLoader::get_ids() const
{
    std::vector<std::string> ids;
    ids.reserve(snp_info_.size());
    for (const auto& s : snp_info_)
    {
        ids.emplace_back(s.id);
    }
    return ids;
}

}  // namespace gelex::detail
