#include "snp_effect_loader.h"

#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "../src/data/parser.h"  // open_file, try_parse_number, parse_string, count_total_lines
#include "gelex/exception.h"

namespace gelex
{

SnpEffects SnpEffectLoader::load(const std::filesystem::path& snp_effect_path)
{
    size_t total_lines = detail::count_total_lines(snp_effect_path);

    auto file = detail::open_file<std::ifstream>(snp_effect_path, std::ios::in);

    std::string line;
    if (!std::getline(file, line))
    {
        throw DataParseException(
            std::format("Empty .snp.eff file: {}", snp_effect_path.string()));
    }

    std::vector<std::string_view> header;
    detail::parse_string(line, header);
    const auto indices = assign_column_indices(header);

    if (!indices.has_required_columns())
    {
        throw HeaderFormatException(
            std::format(
                "Missing required columns (ID, A1, A2, A1Frq, Add) in file: {}",
                snp_effect_path.string()));
    }

    SnpEffects effects;
    if (total_lines > 1)
    {
        effects.reserve(total_lines - 1);
    }

    Eigen::Index idx = 0;
    int line_number = 1;
    const int min_cols_needed = indices.max_required_index() + 1;
    std::vector<std::string_view> row;
    while (std::getline(file, line))
    {
        line_number++;
        if (line.empty())
        {
            continue;
        }

        detail::parse_string(line, row);

        if (static_cast<int>(row.size()) < min_cols_needed)
        {
            throw InconsistentColumnCountException(
                std::format(
                    "Line {} has insufficient columns. Expected at least {}, "
                    "got {}",
                    line_number,
                    min_cols_needed,
                    row.size()));
        }

        try
        {
            auto a1_freq = detail::parse_number<double>(row[indices.a1frq]);
            auto add_val = detail::parse_number<double>(row[indices.add]);

            double dom_val = std::numeric_limits<double>::quiet_NaN();
            if (indices.dom != -1 && indices.dom < static_cast<int>(row.size()))
            {
                dom_val = detail::parse_number<double>(row[indices.dom]);
            }

            SnpEffect effect{
                .index = idx++,
                .A1freq = a1_freq,
                .A1 = row[indices.a1].empty() ? '?' : row[indices.a1][0],
                .A2 = row[indices.a2].empty() ? '?' : row[indices.a2][0],
                .add = add_val,
                .dom = dom_val};

            effects.emplace(std::string(row[indices.id]), effect);
        }
        catch (const GelexException& e)
        {
            throw DataParseException(e.what());
        }
        catch (const std::exception& e)
        {
            throw DataParseException(
                std::format(
                    "Error parsing line {}: {}", line_number, e.what()));
        }
    }

    return effects;
}

bool SnpEffectLoader::has_dom_effects(
    const std::filesystem::path& snp_effect_path)
{
    auto file = detail::open_file<std::ifstream>(snp_effect_path, std::ios::in);

    std::string line;
    if (!std::getline(file, line))
    {
        throw DataParseException("Empty .snp.eff file");
    }

    std::vector<std::string_view> header;
    detail::parse_string(line, header);
    const auto indices = assign_column_indices(header);

    return indices.dom != -1;
}

ColumnIndices SnpEffectLoader::assign_column_indices(
    std::span<const std::string_view> header_columns)
{
    ColumnIndices indices;

    for (int i = 0; i < static_cast<int>(header_columns.size()); ++i)
    {
        const auto& column = header_columns[i];
        if (column == "ID")
        {
            indices.id = i;
        }
        else if (column == "A1")
        {
            indices.a1 = i;
        }
        else if (column == "A2")
        {
            indices.a2 = i;
        }
        else if (column == "A1Frq")
        {
            indices.a1frq = i;
        }
        else if (column == "Add")
        {
            indices.add = i;
        }
        else if (column == "Dom")
        {
            indices.dom = i;
        }
    }

    return indices;
}

}  // namespace gelex
