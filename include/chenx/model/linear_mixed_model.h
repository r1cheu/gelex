#pragma once
#include <armadillo>
#include <vector>

namespace chenx
{
using namespace arma;
class LinearMixedModel
{
   public:
    LinearMixedModel(
        dvec&& y,
        dmat&& X,
        const uvec& z_index,
        dcube&& rands,
        std::vector<std::string>&& rand_names);
    const dvec& y() const { return y_; }
    double y_var() const { return y_var_; }
    const dmat& X() const { return X_; }
    const dvec& beta() const { return beta_; }
    const std::vector<sp_dmat>& Z() const { return Z_; }
    const dcube& zkztr() const { return zkztr_; }
    const std::vector<std::string>& rand_names() const { return rand_names_; }
    const dcube& rands() const { return rands_; }
    const dvec& sigma() const { return sigma_; }

    const double logdet_v() const { return logdet_v_; }
    const dvec& proj_y() const { return proj_y_; }
    const dcube& pdv() const { return pdv_; }
    const dmat& txvx() const { return txvx_; }
    void set_sigma(dvec&& sigma)
    {
        sigma_ = std::move(sigma);
        ComputeV();
        ComputeProj();
        ComputePdV();
    }
    double ComputeLogLikelihood() const;

   private:
    dvec y_;
    double y_var_{};

    dmat X_;
    dvec beta_;

    std::vector<sp_dmat> Z_;
    dcube zkztr_;

    std::vector<std::string> rand_names_;
    dcube rands_;
    dvec sigma_;

    double logdet_v_{};
    dvec proj_y_;
    dmat v_, proj_, txvx_;
    dcube pdv_;

    std::vector<sp_dmat>
    CreateZ(const uword& n_z, const uvec& z_index, const uword& n);
    std::vector<sp_dmat> CreateZ(const uword& n_z, const uword& n);
    dmat ComputeZKZ(const sp_dmat& z, const dmat& k);
    dcube ComputeZKZtR(const std::vector<sp_dmat>& z, const dcube& k);

    void ComputeV();
    void ComputeProj();
    void ComputePdV();
    double VinvLogdet(dmat& V);
};
}  // namespace chenx
