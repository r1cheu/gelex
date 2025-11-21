#pragma once

#include <expected>
#include <filesystem>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "gelex/error.h"

namespace gelex::detail
{

class CovarPredictLoader
{
   public:
    static auto create(const std::filesystem::path& path, bool iid_only)
        -> std::expected<CovarPredictLoader, Error>;

    auto load(const std::unordered_map<std::string, Eigen::Index>& id_map) const
        -> std::map<std::string, std::vector<std::string>>;

    const std::vector<std::string>& names() const { return names_; }
    const std::map<std::string, std::vector<std::string>>& data() const
    {
        return data_;
    }

   private:
    explicit CovarPredictLoader(
        std::vector<std::string>&& names,
        std::map<std::string, std::vector<std::string>>&& data)
        : names_(std::move(names)), data_(std::move(data))
    {
    }

    static auto get_header(std::ifstream& file)
        -> std::expected<std::vector<std::string>, Error>;

    static auto
    read(std::ifstream& file, size_t expected_columns, bool iid_only) -> std::
        expected<std::map<std::string, std::vector<std::string>>, Error>;

    std::vector<std::string> names_;
    std::map<std::string, std::vector<std::string>> data_;
};
}  // namespace gelex::detail
