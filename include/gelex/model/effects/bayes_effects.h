#pragma once
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

#include <armadillo>
#include "base.h"

namespace gelex
{
namespace bayes
{
using arma::uvec;
struct SigmaParam
{
    double nu{-2};
    double s2{0};
    SigmaParam() = default;
};

struct BayesParam
{
    BayesAlphabet type;
    dvec pi;
    uvec pi_num;
    explicit BayesParam(BayesAlphabet type, dvec&& pi, uvec&& pi_num)
        : type(type), pi(std::move(pi)), pi_num(std::move(pi_num))
    {
    }
};

struct BaseEffect
{
    dmat design_mat;
    dvec cols_norm;
    dvec coeff;
    BaseEffect() = default;
    explicit BaseEffect(dmat&& design_matrix);
};

struct FixedEffect : BaseEffect
{
    std::vector<std::string> names;
    std::vector<std::string> levels;
    bool exist{false};

    FixedEffect() = default;
    FixedEffect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        dmat&& design_mat);
};

struct RandomEffect : BaseEffect
{
    std::string name;
    dvec sigma;
    SigmaParam prior;
    RandomEffect(std::string&& name, dvec&& sigma, dmat&& design_matrix);
};

struct GeneticEffect : RandomEffect
{
    BayesParam bayes;
    dvec cols_var;
    dvec u;
    size_t n_zero_var_snp;

    GeneticEffect(
        std::string&& name,
        dvec&& sigma,
        dmat&& design_matrix,
        BayesParam&& bayes);
};

struct Residual
{
    std::string name{"e"};
    double value{0.0};
    SigmaParam prior;
};

struct Mu
{
    std::string name{"mu"};
    double value{0.0};
};

template <typename Effect>
class EffectManager
{
   public:
    void add(Effect&& effect)
    {
        auto index = effects_.size();
        effects_.push_back(std::move(effect));
        index_map_[effects_.back().name] = index;
    }
    Effect* get(const std::string& name)
    {
        auto item = index_map_.find(name);
        return item != index_map_.end() ? &effects_[item->second] : nullptr;
    }
    const Effect* get(const std::string& name) const
    {
        auto item = index_map_.find(name);
        return item != index_map_.end() ? &effects_[item->second] : nullptr;
    }

    size_t size() const { return effects_.size(); }
    bool has_effects() const { return size() > 0; }

    const std::vector<Effect>& effects() const { return effects_; }
    std::vector<std::string> names() const
    {
        std::vector<std::string> names;
        names.clear();
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

    auto begin() { return effects_.begin(); }
    auto end() { return effects_.end(); }
    auto begin() const { return effects_.begin(); }
    auto end() const { return effects_.end(); }

    const Effect& operator[](size_t idx) const { return effects_[idx]; }
    Effect& operator[](size_t idx) { return effects_[idx]; }

   private:
    std::vector<Effect> effects_;
    std::unordered_map<std::string, size_t> index_map_;
};

using RandomEffectManager = EffectManager<RandomEffect>;
using GeneticEffectManager = EffectManager<GeneticEffect>;

template <typename Mat>
dvec sum_square(const Mat& mat)
{
    dvec result(mat.n_cols);

#pragma omp parallel for default(none) shared(result, mat)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        if constexpr (arma::is_Mat_only<Mat>::value)
        {
            result.at(i) = arma::dot(mat.unsafe_col(i), mat.unsafe_col(i));
        }
        else
        {
            result.at(i) = arma::dot(mat.col(i), mat.col(i));
        }
    }
    return result;
}

dvec compute_cols_var(const dmat& mat);

}  // namespace bayes
}  // namespace gelex
