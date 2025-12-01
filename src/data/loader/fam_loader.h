#ifndef GELEX_DATA_LOADER_SAMPLE_LOADER_H
#define GELEX_DATA_LOADER_SAMPLE_LOADER_H

#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

namespace gelex::detail
{

// -----------------------------------------------------------------------------
// FamLoader (Sample Info)
// -----------------------------------------------------------------------------
class FamLoader
{
   public:
    explicit FamLoader(const std::filesystem::path& path, bool iid_only);

    const std::vector<std::string>& ids() const { return ids_; }
    const std::unordered_map<std::string, Eigen::Index>& data() const
    {
        return data_;
    }
    std::vector<std::string>&& take_ids() && { return std::move(ids_); }

   private:
    void set_ids(const std::filesystem::path& path, bool iid_only);
    void set_index_map();
    std::vector<std::string> ids_;
    std::unordered_map<std::string, Eigen::Index> data_;
};

}  // namespace gelex::detail

#endif  // GELEX_DATA_LOADER_SAMPLE_LOADER_H
