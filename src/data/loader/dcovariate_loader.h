#ifndef GELEX_DATA_LOADER_DCOVARIATE_LOADER_H
#define GELEX_DATA_LOADER_DCOVARIATE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>
#include "../src/types/covariates.h"
#include "Eigen/Core"

namespace gelex::detail
{

class DiscreteCovariateLoader
{
   public:
    DiscreteCovariateLoader(const std::filesystem::path& path, bool iid_only);

    [[nodiscard]] DiscreteCovariate load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<std::string>>& data()
        const
    {
        return raw_data_;
    }

   private:
    struct EncodedCovariate
    {
        Eigen::Index active_index = -1;
        Eigen::Index start_col_offset = 0;
    };

    void set_names(std::ifstream& file);
    void set_data(std::ifstream& file, bool iid_only);

    struct IntersectResult
    {
        std::vector<std::string_view> valid_ids;
        std::vector<std::unordered_set<std::string_view>> levels_per_col;
    };
    IntersectResult get_valid_samples_and_levels(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    struct EncodingResult
    {
        std::vector<std::unordered_map<std::string_view, EncodedCovariate>>
            encodings;
        Eigen::Index total_cols;
    };
    EncodingResult build_local_encodings(
        const std::vector<std::unordered_set<std::string_view>>& levels_per_col)
        const;

    DiscreteCovariate build_result(
        const std::unordered_map<std::string, Eigen::Index>& id_map,
        const std::vector<std::string_view>& valid_ids,
        const std::vector<std::unordered_set<std::string_view>>& levels_per_col,
        const EncodingResult& encoding_result) const;

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> raw_data_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_DCOVARIATE_LOADER_H
