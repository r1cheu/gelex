#include "gelex/model/effects.h"

namespace gelex
{

void RandomEffectManager::add_effect(
    std::string&& name,
    effect_type type,
    MatVariant&& design_mat,
    MatVariant&& cov_mat)
{
    auto index = effects_.size();
    effects_.emplace_back(
        std::move(name), type, std::move(design_mat), std::move(cov_mat), 0, 0);
    index_map_[effects_.back().name] = index;

    switch (type)
    {
        case effect_type::random:
            n_random_effects_++;
            random_indices_.emplace_back(index);
            break;
        case effect_type::genetic:
            n_genetic_effects_++;
            genetic_indices_.emplace_back(index);
            break;
        case effect_type::gxe:
            n_gxe_effects_++;
            gxe_indices_.emplace_back(index);
            break;
        case effect_type::residual:
            break;
    }
}

void RandomEffectManager::clear()
{
    effects_.clear();
    index_map_.clear();
    n_random_effects_ = 0;
    n_genetic_effects_ = 0;
    n_gxe_effects_ = 0;
    genetic_indices_.clear();
    random_indices_.clear();
    gxe_indices_.clear();
}
}  // namespace gelex
