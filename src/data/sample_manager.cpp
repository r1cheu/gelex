#include "gelex/data/sample_manager.h"

#include <algorithm>
#include <filesystem>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

#include "gelex/logger.h"
#include "loader/fam_loader.h"

namespace gelex
{

SampleManager::SampleManager(
    const std::filesystem::path& fam_path,
    bool iid_only)
{
    detail::FamLoader fam_loader(fam_path, iid_only);

    common_ids_ = std::move(fam_loader).take_ids();

    std::ranges::sort(common_ids_);

    const auto [first, last] = std::ranges::unique(common_ids_);
    common_ids_.erase(first, last);
}

void SampleManager::intersect(std::span<const std::string_view> ids)
{
    if (common_ids_.empty())
    {
        return;
    }
    if (ids.empty())
    {
        common_ids_.clear();
        return;
    }

    std::vector<std::string_view> sorted_input(ids.begin(), ids.end());
    std::ranges::sort(sorted_input);

    std::vector<std::string> result;
    result.reserve(std::min(common_ids_.size(), sorted_input.size()));

    std::set_intersection(
        common_ids_.begin(),
        common_ids_.end(),
        sorted_input.begin(),
        sorted_input.end(),
        std::back_inserter(result),
        std::less<>{});

    common_ids_ = std::move(result);
}

void SampleManager::finalize()
{
    auto logger = gelex::logging::get();
    logger->info(
        "{} common samples available for analysis after intersection.",
        common_ids_.size());

    common_id_map_.clear();
    common_id_map_.reserve(common_ids_.size());

    for (auto&& [idx, id] : common_ids_ | std::views::enumerate)
    {
        common_id_map_.emplace(id, static_cast<Eigen::Index>(idx));
    }
}

}  // namespace gelex
