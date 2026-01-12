#ifndef GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
#define GELEX_DATA_LOADER_QCOVARIATE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include "../src/types/covariates.h"
#include "Eigen/Core"

namespace gelex::detail
{

class QuantitativeCovariateLoader
{
   public:
    QuantitativeCovariateLoader(
        const std::filesystem::path& path,
        bool iid_only);

    [[nodiscard]] QuantitativeCovariate load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<double>>& data() const
    {
        return data_;
    }

   private:
    void set_names(std::ifstream& file);
    void set_data(std::ifstream& file, bool iid_only);
    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<double>> data_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_QCOVARIATE_LOADER_H
