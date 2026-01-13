#ifndef GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
#define GELEX_DATA_LOADER_QCOVARIATE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

#include "../src/types/covariates.h"

namespace gelex::detail
{

class QuantitativeCovariateLoader
{
   public:
    QuantitativeCovariateLoader(
        const std::filesystem::path& path,
        bool iid_only);

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> QuantitativeCovariate;

    auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

    auto column_names() const -> const std::vector<std::string>&
    {
        return column_names_;
    }

   private:
    static constexpr size_t IdColumnCount = 2;

    std::vector<std::string> column_names_;
    std::vector<std::string> sample_ids_;
    std::vector<std::vector<double>> columns_;

    auto init_columns(std::ifstream& file) -> void;
    auto fill_columns(std::ifstream& file, bool iid_only) -> void;
    [[nodiscard]] static auto is_valid_covariate_value(double value) -> bool;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
