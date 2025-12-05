#include "snp_effect_loader.h"

#include <format>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include <Eigen/Core>

#include "../src/data/parser.h"  // open_file, try_parse_number, parse_string, count_total_lines
#include "gelex/exception.h"

namespace gelex::detail
{

SnpEffectLoader::SnpEffectLoader(const std::filesystem::path& snp_effect_path)
{
    try
    {
        load(snp_effect_path);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", snp_effect_path.string(), e.what()));
    }
}

void SnpEffectLoader::load(const std::filesystem::path& snp_effect_path)
{
    auto file = detail::open_file<std::ifstream>(snp_effect_path, std::ios::in);

    std::string line;
    std::getline(file, line);
    ColumnIndices indices;
    parse_header(line, indices);

    if (!indices.has_required_columns())
    {
        throw HeaderFormatException(
            std::format("missing required columns (ID, A1, A2, A1Frq, Add)"));
    }

    snp_effects_.reserve(1024);
    current_index_ = 0;
    has_dom_ = (indices.dom != -1);

    int line_number = 1;
    while (std::getline(file, line))
    {
        line_number++;
        if (line.empty())
        {
            continue;
        }
        parse_line(line, line_number, indices);
    }
}

void SnpEffectLoader::parse_header(
    std::string_view line,
    ColumnIndices& indices)
{
    std::vector<std::string_view> header;
    detail::parse_string(line, header);
    indices = assign_column_indices(header);
}

void SnpEffectLoader::parse_line(
    std::string_view line,
    int line_number,
    const ColumnIndices& indices)
{
    std::vector<std::string_view> row;
    detail::parse_string(line, row);

    const int min_cols_needed = indices.max_required_index() + 1;
    if (static_cast<int>(row.size()) < min_cols_needed)
    {
        throw InconsistentColumnCountException(
            std::format(
                "{}: has insufficient columns. Expected at least {}, got "
                "{}",
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
            .index = current_index_++,
            .A1freq = a1_freq,
            .A1 = row[indices.a1].empty() ? '?' : row[indices.a1][0],
            .A2 = row[indices.a2].empty() ? '?' : row[indices.a2][0],
            .add = add_val,
            .dom = dom_val};

        snp_effects_.emplace(std::string(row[indices.id]), effect);
    }
    catch (const GelexException& e)
    {
        throw DataParseException(std::format("{}: {}", line_number, e.what()));
    }
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

bool has_dom_effect_column(const std::filesystem::path& snp_effect_path)
{
    auto file = detail::open_file<std::ifstream>(snp_effect_path, std::ios::in);
    std::string line;
    if (!std::getline(file, line))
    {
        throw FileFormatException(
            std::format("{}: empty file", snp_effect_path.string()));
    }

    std::vector<std::string_view> header;
    detail::parse_string(line, header);

    for (const auto& column : header)
    {
        if (column == "Dom")
        {
            return true;
        }
    }
    return false;
}

}  // namespace gelex::detail
