#ifndef GELEX_PREDICTOR_PREDICT_DCOVARIATE_LOADER_H
#define GELEX_PREDICTOR_PREDICT_DCOVARIATE_LOADER_H

#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

class DcovarPredictLoader
{
   public:
    explicit DcovarPredictLoader(
        const std::filesystem::path& path,
        bool iid_only);

    auto load(const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> std::map<std::string, std::vector<std::string>>;

    const std::vector<std::string>& names() const { return names_; }
    const std::unordered_map<std::string, std::vector<std::string>>& data()
        const
    {
        return data_;
    }

   private:
    void set_names(std::ifstream& file);
    void set_data(std::ifstream& file, bool iid_only);

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> data_;
};
}  // namespace gelex::detail

#endif  // GELEX_PREDICTOR_PREDICT_DCOVARIATE_LOADER_H
