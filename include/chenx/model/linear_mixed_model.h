#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include <armadillo>

namespace chenx
{
using arma::dcube;
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
using arma::uvec;
using arma::uword;

struct LinearMixedModelParams
{
    dvec beta;
    dvec sigma;
    std::vector<std::string> individuals;
    std::vector<std::string> dropped_individuals;
};

class LinearMixedModel
{
   public:
    LinearMixedModel(
        dmat&& y,
        dmat&& X,
        dcube&& covar_matrices_rand,
        std::vector<std::string>&& random_effect_names);

    uint64_t num_random_effects() const { return num_random_effects_; }
    uint64_t num_individuals() const { return num_individuals_; }
    uint64_t num_fixed_effects() const { return num_fixed_effects_; }

    const dmat& y() const { return y_; }
    const dmat& X() const { return X_; }
    double y_var() const { return y_var_; }
    const dmat& U() const { return U_; }
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }
    const dvec& proj_y() const { return proj_y_; }
    const dcube& pdv() const { return pdv_; }
    const dmat& v() const { return v_; }
    const dmat& tx_vinv_x() const { return tx_vinv_x_; }
    const dcube& zkzt() const { return zkzt_; }
    const std::vector<std::string>& random_effect_names() const
    {
        return random_effect_names_;
    }
    void set_sigma(dvec&& sigma);
    void set_beta(dvec&& beta) { beta_ = std::move(beta); }
    void set_U(dmat&& U) { U_ = std::move(U); }

    double ComputeLogLikelihood() const;
    void Reset();

   private:
    uint64_t num_random_effects_{};
    uint64_t num_individuals_{};
    uint64_t num_fixed_effects_{};

    dmat y_;
    double y_var_{};

    dmat X_;
    dmat U_;
    dvec beta_;

    dcube zkzt_;

    std::vector<std::string> random_effect_names_;
    dvec sigma_;

    double logdet_v_{};
    dvec proj_y_;
    dmat v_, proj_, tx_vinv_x_;
    dcube pdv_;

    void ComputeV();
    void ComputeProj();
    void ComputePdV();
    static double VinvLogdet(dmat& V);
};
}  // namespace chenx
