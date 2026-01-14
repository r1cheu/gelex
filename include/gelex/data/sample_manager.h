#pragma once

#include <filesystem>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

class SampleManager
{
   public:
    explicit SampleManager(
        const std::filesystem::path& fam_path,
        bool iid_only = false);

    SampleManager(const SampleManager&) = delete;
    SampleManager(SampleManager&&) noexcept = default;
    SampleManager& operator=(const SampleManager&) = delete;
    SampleManager& operator=(SampleManager&&) noexcept = default;
    ~SampleManager() = default;

    void intersect(std::span<const std::string> ids);

    void finalize();

    [[nodiscard]] const std::vector<std::string>& common_ids() const
    {
        return common_ids_;
    }

    [[nodiscard]] const std::unordered_map<std::string, Eigen::Index>&
    common_id_map() const
    {
        return common_id_map_;
    }

    [[nodiscard]] size_t num_common_samples() const
    {
        return common_ids_.size();
    }
    [[nodiscard]] bool has_common_samples() const
    {
        return !common_ids_.empty();
    }

   private:
    std::vector<std::string> common_ids_;
    std::unordered_map<std::string, Eigen::Index> common_id_map_;
};

}  // namespace gelex
