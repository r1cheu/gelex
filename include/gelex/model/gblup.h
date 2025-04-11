#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include <armadillo>

namespace gelex
{
using arma::dcube;
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
using arma::uvec;
using arma::uword;

class GBLUP
{
   public:
    GBLUP(
        dmat y,
        dmat X,
        dcube covar_matrices_rand,
        std::vector<std::string> random_effect_names);

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
    void set_sigma(dvec sigma);
    void set_beta(dvec beta) { beta_ = std::move(beta); }
    void set_U(dmat U) { U_ = std::move(U); }

    double computeLogLikelihood() const;
    void reset();

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
    std::vector<sp_dmat> group_eff;
    dvec sigma_;

    double logdet_v_{};
    dvec proj_y_;
    dmat v_, proj_, tx_vinv_x_;
    dcube pdv_;

    void computeV();
    void computeProj();
    void computePdV();
    static double VinvLogdet(dmat& V);
};

class GBLUPParams
{
   public:
    GBLUPParams(
        dvec beta,
        dvec sigma,
        dvec proj_y,
        std::vector<std::string> dropped_individuals);
    GBLUPParams(
        const GBLUP& model,
        std::vector<std::string> dropped_individuals);
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }
    const dvec& proj_y() const { return proj_y_; }
    const std::vector<std::string>& dropped_individuals() const
    {
        return dropped_individuals_;
    }

    void set_beta(dvec beta) { beta_ = std::move(beta); }
    void set_sigma(dvec sigma) { sigma_ = std::move(sigma); }
    void set_proj_y(dvec proj_y) { proj_y_ = std::move(proj_y); }
    void set_dropped_individuals(std::vector<std::string> dropped_individuals)
    {
        dropped_individuals_ = std::move(dropped_individuals);
    }

   private:
    dvec beta_;
    dvec sigma_;
    dmat X_;
    dvec proj_y_;
    std::vector<std::string> dropped_individuals_;
};

}  // namespace gelex
