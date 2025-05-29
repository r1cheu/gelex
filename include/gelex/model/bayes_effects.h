#pragma once
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include <armadillo>

#include "base_effects.h"

namespace gelex
{
using arma::uvec;

struct SigmaPrior
{
    double nu{-2};
    double s2{0};
};

struct Pi
{
    dvec prop;
    uvec count;
};

struct BaseEffectDesign
{
    explicit BaseEffectDesign(dmat&& design_mat_);

    dmat design_mat;
    dvec cols_norm;
};

struct BaseEffectState
{
    explicit BaseEffectState(size_t n_coeff);
    dvec coeff;
};

struct FixedEffectDesign : BaseEffectDesign
{
    FixedEffectDesign(
        std::vector<std::string>&& names_,
        std::vector<std::string>&& levels_,
        dmat&& design_mat_);

    std::vector<std::string> names;
    std::vector<std::string> levels;
};

struct FixedEffectState : BaseEffectState
{
    using BaseEffectState::BaseEffectState;
    explicit operator bool() const { return !coeff.empty(); }
};

struct RandomEffectDesign : BaseEffectDesign
{
    RandomEffectDesign(std::string&& name_, dmat&& design_mat_);

    std::string name;
    SigmaPrior prior;
};

struct RandomEffectState : BaseEffectState
{
    explicit RandomEffectState(size_t n_coeff);
    dvec coeff;
    dvec sigma{0};  // set to dvec (not double) for consistency
};

struct GeneticEffectDesign : RandomEffectDesign
{
    GeneticEffectDesign(
        std::string&& name_,
        dmat&& design_mat_,
        BayesAlphabet type_,
        dvec&& sigma_,
        dvec&& pi_);

    dvec cols_var;
    size_t n_zero_var_snp;
    BayesAlphabet type;
    dvec pi;
    dvec sigma;
};

struct GeneticEffectState : BaseEffectState
{
    explicit GeneticEffectState(
        size_t n_individual,
        size_t n_coeff,
        const dvec& pi_prop,
        const dvec& sigma_);
    dvec u;
    Pi pi;
    dvec sigma;
};

struct Residual
{
    std::string name{"e"};
    SigmaPrior prior;
    double value{0.0};
};

struct Mu
{
    std::string name{"mu"};
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
dvec sum_square(const Mat& mat)
{
    dvec result(mat.n_cols);

#pragma omp parallel for default(none) shared(result, mat)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        result.at(i) = arma::dot(mat.col(i), mat.col(i));
    }
    return result;
}

dvec compute_cols_var(const dmat& mat);

}  // namespace gelex
