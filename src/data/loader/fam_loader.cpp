#include "fam_loader.h"

#include <filesystem>
#include <fstream>
#include <ranges>
#include <vector>

#include <Eigen/Core>

#include <fmt/ranges.h>
#include "../parser.h"
#include "gelex/exception.h"

namespace gelex::detail
{

FamLoader::FamLoader(const std::filesystem::path& path, bool iid_only)
{
    try
    {
        set_ids(path, iid_only);
        set_index_map();
    }
    catch (const GelexException& e)
    {
        throw FileFormatException(
            std::format("{}:{}", path.string(), e.what()));
    }
}

void FamLoader::set_ids(const std::filesystem::path& path, bool iid_only)
{
    ids_.clear();
    auto file = detail::open_file<std::ifstream>(path, std::ios::in);
    char delimiter = detect_file_delimiter(file);

    ids_.reserve(1024);  // Start small
    std::string line;
    int n_line = 0;

    while (std::getline(file, line))
    {
        n_line++;
        if (line.empty())
        {
            continue;
        }
        try
        {
            ids_.emplace_back(parse_id(line, iid_only, delimiter));
        }
        catch (const GelexException& e)
        {
            throw DataParseException(std::format("{}: {}", n_line, e.what()));
        }
    }
}

void FamLoader::set_index_map()
{
    data_.reserve(ids_.size());
    for (auto&& [idx, id] : ids_ | std::views::enumerate)
    {
        data_[id] = static_cast<Eigen::Index>(idx);
    }
}

}  // namespace gelex::detail
