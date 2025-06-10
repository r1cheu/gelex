#pragma once
#include <cstddef>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include <armadillo>

#include "base_effects.h"

namespace gelex
{
using MatVariant = std::variant<arma::dmat, arma::sp_dmat>;

struct FixedEffect
{
    std::vector<std::string> names;
    std::vector<std::string> levels;
    MatVariant design_mat;
    arma::dvec beta;

    size_t size() const { return beta.n_elem; }
    void clear();
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
    void add(
        std::string&& name,
        effect_type type,
        MatVariant&& design_mat,
        MatVariant&& cov_mat);

    RandomEffect* get(const std::string& name);
    const RandomEffect* get(const std::string& name) const;

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

    void set_sigma(const arma::dvec& sigma);

    void set_se(const arma::dvec& se);

    const arma::dvec& sigma() const { return sigma_; }
    void set_hess_inv(const arma::dmat& hess_inv) { hess_inv_ = hess_inv; }
    const arma::dmat& hess_inv() const { return hess_inv_; }

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

    arma::dmat hess_inv_;
    arma::dvec sigma_;
};
}  // namespace gelex
