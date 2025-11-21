#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include "gelex/error.h"

namespace gelex
{

class SampleManager
{
   public:
    static auto create(
        const std::filesystem::path& fam_path,
        bool iid_only = false) -> std::expected<SampleManager, Error>;

    SampleManager(const SampleManager&) = delete;
    SampleManager(SampleManager&&) noexcept = default;
    SampleManager& operator=(const SampleManager&) = delete;
    SampleManager& operator=(SampleManager&&) noexcept = default;
    ~SampleManager() = default;

    void intersect(std::span<const std::string_view> ids);
    void finalize();

    const std::vector<std::string>& common_ids() const { return common_ids_; }
    const std::unordered_map<std::string, Eigen::Index>& common_id_map() const
    {
        return common_id_map_;
    }
    size_t num_common_samples() const { return common_ids_.size(); }
    bool has_common_samples() const { return !common_ids_.empty(); }

   private:
    SampleManager() = default;

    std::vector<std::string> common_ids_;
    std::unordered_set<std::string> common_ids_set_;
    std::unordered_map<std::string, Eigen::Index> common_id_map_;
};

}  // namespace gelex
