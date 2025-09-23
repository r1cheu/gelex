#pragma once

#include <expected>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
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

    // Genotyped sample information from FAM file (for BED file reading)
    const std::vector<std::string>& genotyped_sample_ids() const
    {
        return fam_loader_->ids();
    }
    const std::unordered_map<std::string, Eigen::Index>& genotyped_sample_map()
        const
    {
        return fam_loader_->data();
    }

    size_t num_common_samples() const { return common_ids_.size(); }
    size_t num_genotyped_samples() const
    {
        return fam_loader_ ? fam_loader_->ids().size() : 0;
    }

    bool has_common_samples() const { return !common_ids_.empty(); }
    bool has_genotyped_samples() const { return fam_loader_ != nullptr; }

   private:
    SampleManager() = default;

    std::unique_ptr<detail::FamLoader> fam_loader_;

    std::vector<std::string> common_ids_;
    std::unordered_set<std::string> common_ids_set_;
    std::unordered_map<std::string, Eigen::Index> common_id_map_;
};

}  // namespace gelex
