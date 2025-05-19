#pragma once
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <armadillo>

namespace gelex
{

using arma::dmat;
using arma::dvec;
using MatVariant = std::variant<arma::dmat, arma::sp_dmat>;

enum class effect_type : uint8_t
{
    random,
    genetic,
    gxe,
    residual,
};

struct FixedEffect
{
    std::vector<std::string> names;
    std::vector<std::string> levels;
    MatVariant design_mat;
    dvec beta;

    size_t size() const { return beta.n_elem; }
};

struct RandomEffect
{
    std::string name;
    effect_type type;
    MatVariant design_mat;
    MatVariant cov_mat;
    double sigma;
    double se;
};

class RandomEffectManager
{
   public:
    void add_effect(
        std::string&& name,
        effect_type type,
        MatVariant&& design_mat,
        MatVariant&& cov_mat);

    RandomEffect* get(const std::string& name)
    {
        auto item = index_map_.find(name);
        return item != index_map_.end() ? &effects_[item->second] : nullptr;
    }
    const RandomEffect* get(const std::string& name) const
    {
        auto item = index_map_.find(name);
        return item != index_map_.end() ? &effects_[item->second] : nullptr;
    }

    size_t size() const { return effects_.size(); }

    size_t n_random_effects() const { return n_random_effects_; }
    size_t n_genetic_effects() const { return n_genetic_effects_; }
    size_t n_gxe_effects() const { return n_gxe_effects_; }

    bool has_random_effects() const { return n_random_effects_ > 0; }
    bool has_genetic_effects() const { return n_genetic_effects_ > 0; }
    bool has_gxe_effects() const { return n_gxe_effects_ > 0; }

    const std::vector<size_t>& genetic_indices() const
    {
        return genetic_indices_;
    }
    const std::vector<size_t>& random_indices() const
    {
        return random_indices_;
    }
    const std::vector<size_t>& gxe_indices() const { return gxe_indices_; }
    size_t residul_index() const { return residual_index_; }

    const std::vector<RandomEffect>& effects() const { return effects_; }

    void clear();

    auto begin() { return effects_.begin(); }
    auto end() { return effects_.end(); }
    auto begin() const { return effects_.begin(); }
    auto end() const { return effects_.end(); }

    const RandomEffect& operator[](size_t idx) const { return effects_[idx]; }
    RandomEffect& operator[](size_t idx) { return effects_[idx]; }

    void set_sigma(const dvec& sigma)
    {
        sigma_ = sigma;
        for (size_t i = 0; i < effects_.size(); ++i)
        {
            effects_[i].sigma = sigma_.at(i);
        }
    }

    void set_se(const dvec& se)
    {
        for (size_t i = 0; i < effects_.size(); ++i)
        {
            effects_[i].se = se.at(i);
        }
    }

    const dvec& sigma() const { return sigma_; }
    void set_hess_inv(const dmat& hess_inv) { hess_inv_ = hess_inv; }
    const dmat& hess_inv() const { return hess_inv_; }

   private:
    std::vector<RandomEffect> effects_;
    std::unordered_map<std::string, size_t> index_map_;

    size_t n_random_effects_{};
    size_t n_genetic_effects_{};
    size_t n_gxe_effects_{};

    std::vector<size_t> random_indices_;
    std::vector<size_t> genetic_indices_;
    std::vector<size_t> gxe_indices_;
    size_t residual_index_{};

    dmat hess_inv_;
    dvec sigma_;
};

}  // namespace gelex
