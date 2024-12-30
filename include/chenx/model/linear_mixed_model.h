#pragma once
#include <memory>
#include <vector>

#include <spdlog/logger.h>
#include <armadillo>

namespace chenx
{
using arma::dcube;
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
using arma::uvec;
using arma::uword;

class LinearMixedModel
{
   public:
    LinearMixedModel(
        dvec&& y,
        dmat&& X,
        std::vector<sp_dmat>&& Z,
        dcube&& covar_matrices_rand,
        std::vector<std::string>&& rand_names);
    const dvec& y() const { return y_; }
    double y_var() const { return y_var_; }
    const dmat& X() const { return X_; }
    const dvec& beta() const { return beta_; }
    const std::vector<sp_dmat>& Z() const { return Z_; }
    const dcube& zkztr() const { return zkztr_; }
    const dmat& v() const { return v_; }
    const std::vector<std::string>& rand_names() const { return rand_names_; }
    const dcube& covar_matrices_rand() const { return covar_matrices_rand_; }
    const dvec& sigma() const { return sigma_; }

    double logdet_v() const { return logdet_v_; }
    const dvec& proj_y() const { return proj_y_; }
    const dcube& pdv() const { return pdv_; }
    const dmat& tx_vinv_x() const { return tx_vinv_x_; }
    void set_sigma(dvec&& sigma)
    {
        sigma_ = std::move(sigma);
        ComputeV();
        ComputeProj();
        ComputePdV();
    }
    void set_beta(dvec&& beta) { beta_ = std::move(beta); }
    double ComputeLogLikelihood() const;

   private:
    dvec y_;
    double y_var_{};

    dmat X_;
    dvec beta_;

    std::vector<sp_dmat> Z_;
    dcube zkztr_;

    std::vector<std::string> rand_names_;
    dcube covar_matrices_rand_;
    dvec sigma_;

    double logdet_v_{};
    dvec proj_y_;
    dmat v_, proj_, tx_vinv_x_;
    dcube pdv_;
    std::shared_ptr<spdlog::logger> logger_;

    static dmat ComputeZKZ(const sp_dmat& z, const dmat& k);
    dcube ComputeZKZtR();

    void ComputeV();
    void ComputeProj();
    void ComputePdV();
    double VinvLogdet(dmat& V);
};
}  // namespace chenx
