#include "gelex/model/effects.h"

namespace gelex
{

void GroupEffectManager::add_effect(
    std::string name,
    effect_type type,
    MatVariant design_mat,
    MatVariant cov_mat)
{
    auto it = index_map_.find(name);
    auto index = effects_.size();
    if (it == index_map_.end())
    {
        effects_.push_back(
            {std::move(name),
             type,
             std::move(design_mat),
             std::move(cov_mat),
             0});
        index_map_[effects_.back().name] = index;
    }
    else
    {
        effects_[it->second].design_mat = std::move(design_mat);
        effects_[it->second].cov_mat = std::move(cov_mat);
    }
    switch (type)
    {
        case effect_type::group:
            n_group_effects_++;
            group_indices_.emplace_back(index);
            break;
        case effect_type::genetic:
            n_genetic_effects_++;
            genetic_indices_.emplace_back(index);
            break;
        case effect_type::gxe:
            n_gxe_effects_++;
            gxe_indices_.emplace_back(index);
            break;
    }
}

void GroupEffectManager::clear()
{
    effects_.clear();
    index_map_.clear();
    n_group_effects_ = 0;
    n_genetic_effects_ = 0;
    n_gxe_effects_ = 0;
    genetic_indices_.clear();
    group_indices_.clear();
    gxe_indices_.clear();
}
}  // namespace gelex
