#pragma once
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

class LinearMixedModel
{
   public:
    LinearMixedModel(
        dmat&& y,
        dmat&& X,
        dcube&& covar_matrices_rand,
        std::vector<std::string>&& rand_names);
    const dmat& y() const { return y_; }
    const dmat& X() const { return X_; }
    double y_var() const { return y_var_; }
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }
    const dvec& proj_y() const { return proj_y_; }
    const dcube& pdv() const { return pdv_; }
    const dmat& v() const { return v_; }
    const dmat& tx_vinv_x() const { return tx_vinv_x_; }
    const std::vector<std::string>& rand_names() const { return rand_names_; }

    void Reset()
    {
        auto n_rand = rand_names_.size();
        set_sigma(dvec(
            n_rand, arma::fill::value(y_var_ / static_cast<double>(n_rand))));
        set_beta(dvec(X_.n_cols, arma::fill::zeros));
    };
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
    dmat y_;
    double y_var_{};

    dmat X_;
    dvec beta_;

    dcube zkzt_;
    sp_dmat r_;

    std::vector<std::string> rand_names_;
    dvec sigma_;

    double logdet_v_{};
    dvec proj_y_;
    dmat v_, proj_, tx_vinv_x_;
    dcube pdv_;

    // static dmat ComputeZKZ(const sp_dmat& z, const dmat& k);
    //  dcube ComputeZKZtR(dcube&& covar_matrices_rand);

    void ComputeV();
    void ComputeProj();
    void ComputePdV();
    double VinvLogdet(dmat& V);
};
}  // namespace chenx
