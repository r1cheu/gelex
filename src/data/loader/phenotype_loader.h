#ifndef GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
#define GELEX_DATA_LOADER_PHENOTYPE_LOADER_H

#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>

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

    [[nodiscard]] Eigen::VectorXd load(
        const std::unordered_map<std::string, Eigen::Index>& id_map) const;

    const std::string& name() const { return name_; }
    const std::unordered_map<std::string, double>& data() const
    {
        return data_;
    }

   private:
    void set_name(std::ifstream& file, int pheno_column);
    void set_data(std::ifstream& file, int pheno_column, bool iid_only);
    std::string name_;
    std::unordered_map<std::string, double> data_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_PHENOTYPE_LOADER_H
