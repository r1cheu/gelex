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

PhenotypeLoader::PhenotypeLoader(
    const std::filesystem::path& path,
    int pheno_column,
    bool iid_only)
{
    auto file = open_file<std::ifstream>(path, std::ios::in);
    try
    {
        init_column(file, pheno_column);
        fill_data(file, pheno_column, iid_only);
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}

auto PhenotypeLoader::init_column(std::ifstream& file, int pheno_column) -> void
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

auto PhenotypeLoader::fill_data(
    std::ifstream& file,
    int pheno_column,
    bool iid_only) -> void
{
    std::string line;
    size_t n_line = 1;

    sample_ids_.reserve(1024);
    values_.reserve(1024);

    while (std::getline(file, line))
    {
        ++n_line;
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
            sample_ids_.emplace_back(parse_id(line, iid_only));
            values_.push_back(value);
        }
        catch (const GelexException& e)
        {
            throw DataParseException(std::format("{}: {}", n_line, e.what()));
        }
    }
}

auto PhenotypeLoader::load(
    const std::unordered_map<std::string, Eigen::Index>& id_map) const
    -> Eigen::VectorXd
{
    Eigen::VectorXd result(id_map.size());
    result.setConstant(std::numeric_limits<double>::quiet_NaN());

    for (size_t i = 0; i < sample_ids_.size(); ++i)
    {
        if (auto it = id_map.find(sample_ids_[i]); it != id_map.end())
        {
            result(it->second) = values_[i];
        }
    }
    return result;
}

}  // namespace gelex::detail
