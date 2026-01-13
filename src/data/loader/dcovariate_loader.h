#ifndef GELEX_DATA_LOADER_DCOVARIATE_LOADER_H
#define GELEX_DATA_LOADER_DCOVARIATE_LOADER_H

#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "../src/types/covariates.h"

namespace gelex::detail
{

class DiscreteCovariateLoader
{
   public:
    DiscreteCovariateLoader(const std::filesystem::path& path, bool iid_only);

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> DiscreteCovariate;

    auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

    auto column_names() const -> const std::vector<std::string>&
    {
        return column_names_;
    }

   private:
    struct StringHash
    {
        using is_transparent = void;
        auto operator()(std::string_view sv) const noexcept -> size_t
        {
            return std::hash<std::string_view>{}(sv);
        }
    };

    struct ColumnData
    {
        std::vector<std::string> levels;
        std::vector<uint16_t> data;
        std::unordered_map<std::string, uint16_t, StringHash, std::equal_to<>>
            level_map;
    };

    std::vector<std::string> column_names_;
    std::vector<std::string> sample_ids_;
    std::vector<ColumnData> columns_;

    // set column names and initialize column data structures
    auto init_columns(std::ifstream& file) -> void;
    auto fill_columns(std::ifstream& file, bool iid_only) -> void;
    static auto get_or_add_level(ColumnData& column, std::string_view level)
        -> uint16_t;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_DCOVARIATE_LOADER_H
