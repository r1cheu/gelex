#include "gelex/data/sample_manager.h"

#include <algorithm>
#include <expected>
#include <filesystem>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "../src/data/loader.h"
#include "gelex/error.h"
#include "gelex/logger.h"

namespace gelex
{

auto SampleManager::create(const std::filesystem::path& fam_path, bool iid_only)
    -> std::expected<SampleManager, Error>
{
    SampleManager manager;

    auto fam_loader = detail::FamLoader::create(fam_path, iid_only);
    if (!fam_loader)
    {
        return std::unexpected(fam_loader.error());
    }

    std::unordered_set<std::string> common_ids_set;
    for (const auto& [id, _] : fam_loader.value().data())
    {
        common_ids_set.insert(id);
    }

    manager.fam_loader_
        = std::make_unique<detail::FamLoader>(std::move(*fam_loader));
    manager.common_ids_set_ = std::move(common_ids_set);

    return manager;
}

void SampleManager::intersect(std::span<const std::string_view> ids)
{
    std::unordered_set<std::string_view> id_set{ids.begin(), ids.end()};
    std::erase_if(
        common_ids_set_, [&](const auto& id) { return !id_set.contains(id); });
}

void SampleManager::finalize()
{
    std::vector<std::string> id_vec{
        std::make_move_iterator(common_ids_set_.begin()),
        std::make_move_iterator(common_ids_set_.end())};
    common_ids_ = std::move(id_vec);
    std::ranges::sort(common_ids_);

    auto logger = gelex::logging::get();
    logger->info(
        "{} common samples available for analysis after intersection.",
        common_ids_.size());

    common_id_map_.clear();
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(common_ids_.size());
         ++i)
    {
        common_id_map_.emplace(common_ids_[i], i);
    }
}

}  // namespace gelex
