#include "gelex/model/effects.h"
#include "gelex/model/freq/effects.h"

namespace gelex
{
using arma::dvec;

void FixedEffect::clear()
{
    auto clear = [](auto& mat) { mat.clear(); };
    names.clear();
    levels.clear();
    std::visit(clear, design_mat);
    beta.clear();
}

RandomEffect* RandomEffectManager::get(const std::string& name)
{
    auto item = index_map_.find(name);
    return item != index_map_.end() ? &effects_[item->second] : nullptr;
}

const RandomEffect* RandomEffectManager::get(const std::string& name) const
{
    auto item = index_map_.find(name);
    return item != index_map_.end() ? &effects_[item->second] : nullptr;
}

void RandomEffectManager::add(
    std::string&& name,
    effect_type type,
    MatVariant&& design_mat,
    MatVariant&& cov_mat)
{
    auto index = effects_.size();
    effects_.emplace_back(
        std::move(name),
        type,
        std::move(design_mat),
        std::move(cov_mat),
        arma::dvec{},
        arma::dvec{},
        0,
        0);
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

void RandomEffectManager::set_sigma(const dvec& sigma)
{
    sigma_ = sigma;
    for (size_t i = 0; i < effects_.size(); ++i)
    {
        effects_[i].sigma = sigma[i];
    }
}

void RandomEffectManager::set_se(const dvec& se)
{
    for (size_t i = 0; i < effects_.size(); ++i)
    {
        effects_[i].se = se[i];
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
    residual_index_ = 0;

    sigma_.clear();
    hess_inv_.clear();
};
}  // namespace gelex
