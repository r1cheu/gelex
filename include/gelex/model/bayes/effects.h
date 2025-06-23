#pragma once
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include <armadillo>

#include "gelex/dist.h"
#include "gelex/model/effects.h"

namespace gelex
{

struct Pi
{
    arma::dvec prop;
    arma::uvec count;
};

struct BaseEffectDesign
{
    explicit BaseEffectDesign(arma::dmat&& design_mat_);

    arma::dmat design_mat;
    arma::dvec cols_norm;
};

struct BaseEffectState
{
    explicit BaseEffectState(size_t n_coeff);
    arma::dvec coeff;
};

struct FixedEffectDesign : BaseEffectDesign
{
    FixedEffectDesign(
        std::vector<std::string>&& names_,
        std::vector<std::string>&& levels_,
        arma::dmat&& design_mat_);

    std::vector<std::string> names;
    std::vector<std::string> levels;
};

struct FixedEffectState : BaseEffectState
{
    using BaseEffectState::BaseEffectState;
};

struct RandomEffectDesign : BaseEffectDesign
{
    RandomEffectDesign(std::string&& name_, arma::dmat&& design_mat_);

    std::string name;
    ScaledInvChiSqParams prior;
    arma::dvec sigma{0};  // set to dvec (not double) for consistency
};

struct RandomEffectState : BaseEffectState
{
    explicit RandomEffectState(size_t n_coeff, const arma::dvec& init_sigma);
    arma::dvec coeff;
    arma::dvec sigma{0};  // set to dvec (not double) for consistency
};

struct GeneticEffectDesign : RandomEffectDesign
{
    GeneticEffectDesign(
        std::string&& name_,
        arma::dmat&& design_mat_,
        BayesAlphabet type_,
        arma::dvec&& sigma_,
        arma::dvec&& pi_);

    BayesAlphabet type;
    arma::dvec pi;
    arma::dvec sigma;
    arma::dvec mean;
    arma::dvec stddev;
};

struct GeneticEffectState : BaseEffectState
{
    explicit GeneticEffectState(
        BayesAlphabet type_,
        size_t n_individual,
        size_t n_coeff,
        const arma::dvec& pi_prop,
        const arma::dvec& sigma_);

    BayesAlphabet type;
    arma::dvec u;
    Pi pi;
    double genetic_var{};
    double heritability{};
    arma::dvec sigma;
};

struct Residual
{
    std::string name{"e"};
    ScaledInvChiSqParams prior;
    double value{0.0};
};

template <typename Design>
class EffectDesignManager
{
   public:
    void add(Design&& design)
    {
        effects_.push_back(std::move(design));
        index_map_[effects_.back().name] = effects_.size() - 1;
    }

    const Design* get(const std::string& name) const
    {
        if (auto it = index_map_.find(name); it != index_map_.end())
        {
            return &effects_[it->second];
        }
        return nullptr;
    }
    Design* get(const std::string& name)
    {
        if (auto it = index_map_.find(name); it != index_map_.end())
        {
            return &effects_[it->second];
        }
        return nullptr;
    }

    size_t size() const { return effects_.size(); }
    const std::vector<Design>& effects() const { return effects_; }
    std::vector<Design>& effects() { return effects_; }

    std::vector<std::string> names() const
    {
        std::vector<std::string> names;
        names.reserve(effects_.size());
        for (const auto& effect : effects_)
        {
            names.push_back(effect.name);
        }
        return names;
    }

    void clear()
    {
        effects_.clear();
        index_map_.clear();
    }
    const Design& operator[](size_t index) const { return effects_[index]; }
    Design& operator[](size_t index) { return effects_[index]; }
    const Design& back() const { return effects_.back(); }
    Design& back() { return effects_.back(); }

    auto begin() { return effects_.begin(); }
    auto end() { return effects_.end(); }
    auto begin() const { return effects_.begin(); }
    auto end() const { return effects_.end(); }

    explicit operator bool() const { return !effects_.empty(); }

   private:
    std::vector<Design> effects_;
    std::unordered_map<std::string, size_t> index_map_;
};

using RandomEffectDesignManager = EffectDesignManager<RandomEffectDesign>;
using GeneticEffectDesignManager = EffectDesignManager<GeneticEffectDesign>;

std::vector<RandomEffectState> create_thread_states(
    const RandomEffectDesignManager& designs);
std::vector<GeneticEffectState> create_thread_states(
    const GeneticEffectDesignManager& designs);

template <typename Mat>
arma::dvec sum_square(const Mat& mat)
{
    arma::dvec result(mat.n_cols);

#pragma omp parallel for default(none) shared(result, mat)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        result.at(i) = arma::dot(mat.col(i), mat.col(i));
    }
    return result;
}

arma::dvec compute_cols_var(const arma::dmat& mat);

}  // namespace gelex
