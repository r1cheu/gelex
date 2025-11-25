#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "gelex/exception.h"

namespace gelex::detail
{

class CovarPredictLoader
{
   public:
    explicit CovarPredictLoader(
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
    static auto get_header(std::ifstream& file) -> std::vector<std::string>;

    static auto
    read(std::ifstream& file, size_t expected_columns, bool iid_only)
        -> std::unordered_map<std::string, std::vector<std::string>>;

    std::vector<std::string> names_;
    std::unordered_map<std::string, std::vector<std::string>> data_;
};
}  // namespace gelex::detail
