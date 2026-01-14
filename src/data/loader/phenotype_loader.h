#ifndef GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
#define GELEX_DATA_LOADER_PHENOTYPE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

class PhenotypeLoader
{
   public:
    PhenotypeLoader(
        const std::filesystem::path& path,
        int pheno_column,
        bool iid_only);

    [[nodiscard]] auto load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> Eigen::VectorXd;

    auto name() const -> const std::string& { return name_; }

    auto sample_ids() const -> const std::vector<std::string>&
    {
        return sample_ids_;
    }

   private:
    std::string name_;
    std::vector<std::string> sample_ids_;
    std::vector<double> values_;

    auto init_column(std::ifstream& file, int pheno_column) -> void;
    auto fill_data(std::ifstream& file, int pheno_column, bool iid_only)
        -> void;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
