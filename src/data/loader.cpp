#include "gelex/data/loader.h"

#include <algorithm>  // For std::ranges::sort, std::ranges::set_intersection
#include <filesystem>
#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
//
#include <Eigen/Dense>

#include "../src/data/parser.h"

namespace gelex::detail
{

void validate_path_or_throw(const std::filesystem::path& path) {}

std::unordered_set<std::string> get_ids_from_file(
    const std::filesystem::path& path,
    bool iid_only,
    bool skip_header)
{
    auto file = detail::open_or_throw<std::ifstream>(path);
    std::unordered_set<std::string> ids;
    std::string line;

    if (skip_header)
    {
        std::getline(file, line);
    }

    while (std::getline(file, line))
    {
        if (line.empty())
        {
            continue;
        }
        if (auto id = parse_id(line, iid_only); id)
        {
            ids.insert(*id);
        }
    }
    return ids;
}

class PhenotypeLoader
{
   public:
    explicit PhenotypeLoader(
        std::string_view path,
        int pheno_column,
        bool iid_only)
        : path_(path), pheno_column_(pheno_column), iid_only_(iid_only)
    {
    }

    std::unordered_set<std::string> get_sample_ids()
    {
        auto file = detail::open_or_throw<std::ifstream>(path_);
        std::unordered_set<std::string> sample_ids;

        std::string line;
        std::getline(file, line);
        const auto header = parse_header(line, path_)
                            | std::ranges::to<std::vector<std::string>>();
        int n_samples{};
        while (std::getline(file, line))
        {
            if (line.empty())
            {
                continue;
            }
            if (auto id_str = parse_id(line, iid_only_); id_str)
            {
                if (auto value = parse_nth_double(line, pheno_column_))
                {
                    sample_ids.emplace(std::move(*id_str));
                }
            }
            n_samples++;
        }
        return sample_ids;
    }

    Eigen::MatrixXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        bool iid_only) const
    {
        auto file = detail::open_or_throw<std::ifstream>(path_);
        std::string line;

        std::getline(file, line);
        const auto header = detail::parse_header(line, path_)
                            | std::ranges::to<std::vector<std::string>>();
        const auto num_cols = static_cast<Eigen::Index>(header.size());

        Eigen::MatrixXd data_matrix(id_map.size(), num_cols);

        while (std::getline(file, line))
        {
            if (line.empty())
            {
                continue;
            }
            if (auto id = parse_id(line, iid_only); id)
            {
                if (auto it = id_map.find(*id); it != id_map.end())
                {
                }
            }
        }
        return data_matrix;
    }

   private:
    std::filesystem::path path_;
    int pheno_column_;
    bool iid_only_;
};

}  // namespace gelex::detail
