#include "phenotype_loader.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <vector>

#include <fmt/ranges.h>
#include <Eigen/Core>

#include "gelex/exception.h"

// Internal
#include "../src/data/parser.h"

namespace gelex::detail
{

// =============================================================================
// PhenotypeLoader
// =============================================================================

PhenotypeLoader::PhenotypeLoader(
    const std::filesystem::path& path,
    int pheno_column,
    bool iid_only)
{
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    try
    {
        set_name(file, pheno_column);
        set_data(file, pheno_column, iid_only);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}
void PhenotypeLoader::set_name(std::ifstream& file, int pheno_column)
{
    std::string line;
    std::getline(file, line);
    auto header = parse_header(line);
    if (pheno_column < 2 || static_cast<size_t>(pheno_column) >= header.size())
    {
        throw ColumnRangeException(
            std::format("Phenotype column {} is out of range", pheno_column));
    }
    name_ = std::string(header[pheno_column]);
}

void PhenotypeLoader::set_data(
    std::ifstream& file,
    int pheno_column,
    bool iid_only)
{
    std::string line;
    int n_line = 0;
    data_.reserve(1024);

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            double value = parse_nth_double(line, pheno_column);
            if (std::isnan(value) || std::isinf(value))
            {
                continue;
            }
            data_.emplace(parse_id(line, iid_only), value);
        }
        catch (const GelexException& e)
        {
            throw DataParseException(
                std::format("{}: {}", n_line + 1, e.what()));
        }
    }
}

Eigen::VectorXd PhenotypeLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
{
    Eigen::VectorXd result(id_map.size());
    result.setConstant(std::numeric_limits<double>::quiet_NaN());

    for (const auto& [id, value] : data_)
    {
        if (auto it = id_map.find(id); it != id_map.end())
        {
            result(it->second) = value;
        }
    }
    return result;
}

}  // namespace gelex::detail
