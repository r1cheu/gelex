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
using MatVariant = std::variant<arma::dmat, arma::sp_dmat>;
using arma::dvec;

enum class effect_type : uint8_t
{
    group,
    genetic,
    gxe
};

struct Effect
{
    std::string name;
    effect_type type;
    MatVariant design_mat;
    MatVariant cov_mat;
    double sigma;
    double se;
};

class EffectManager
{
   public:
    void add_effect(
        std::string name,
        effect_type type,
        MatVariant design_mat,
        MatVariant cov_mat);

    Effect* get(const std::string& name)
    {
        auto it = index_map_.find(name);
        return it != index_map_.end() ? &effects_[it->second] : nullptr;
    }

    uint64_t size() const
    {
        return effects_.size() + 1;
    }  // add one for residual

    uint64_t n_group_effects() const { return n_group_effects_; }
    uint64_t n_genetic_effects() const { return n_genetic_effects_; }
    uint64_t n_gxe_effects() const { return n_gxe_effects_; }
    bool has_group_effects() const { return n_group_effects_ > 0; }
    bool has_genetic_effects() const { return n_genetic_effects_ > 0; }
    bool has_gxe_effects() const { return n_gxe_effects_ > 0; }

    const std::vector<uint64_t>& genetic_indices() const
    {
        return genetic_indices_;
    }
    const std::vector<uint64_t>& group_indices() const
    {
        return group_indices_;
    }
    const std::vector<uint64_t>& gxe_indices() const { return gxe_indices_; }
    const std::vector<Effect>& effects() const { return effects_; }
    const dmat& effect_cov() const { return effect_cov_; }
    void set_effect_cov(dmat&& cov)
    {
        effect_cov_ = std::move(cov);
        compute_se();
    }

    void set_residual(double res) { residual_ = res; }
    double residual() const { return residual_; }
    double residual_se() const { return residual_se_; }
    void clear();
    auto begin() { return effects_.begin(); }
    auto end() { return effects_.end(); }
    auto begin() const { return effects_.begin(); }
    auto end() const { return effects_.end(); }

    const Effect& operator[](size_t i) const { return effects_[i]; }
    Effect& operator[](size_t i) { return effects_[i]; }

   private:
    void compute_se()
    {
        dvec se = arma::sqrt(arma::diagvec(-effect_cov_));
        for (size_t i = 0; i < effects_.size(); ++i)
        {
            effects_[i].se = se[i];
        }
        residual_se_ = se.back();
    }
    std::vector<Effect> effects_;
    std::unordered_map<std::string, size_t> index_map_;
    uint64_t n_group_effects_{};
    uint64_t n_genetic_effects_{};
    uint64_t n_gxe_effects_{};

    double residual_{};
    double residual_se_{};

    dmat effect_cov_;

    std::vector<uint64_t> genetic_indices_;
    std::vector<uint64_t> group_indices_;
    std::vector<uint64_t> gxe_indices_;
};

}  // namespace gelex
