#pragma once
#include <armadillo>

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;
constexpr static double add_h2 = 0.5;

struct sigma_prior
{
    double nu;
    double s2;
};

class SimgaPriors
{
   public:
    SimgaPriors() = default;
    explicit SimgaPriors(const dvec& pi) : pi_{pi} {}

    sigma_prior& sigma_g() { return sigma_g_; }
    sigma_prior& sigma_r() { return sigma_r_; }
    sigma_prior& sigma_e() { return sigma_e_; }

    void set_genomic_effect_s2(const dvec& y, const dmat& genotype_mat);

    dvec pi() const { return pi_; }
    double add_h2() const { return add_h2_; }

    void set_pi(dvec pi) { pi_ = std::move(pi); }
    void set_add_h2(double add_h2) { add_h2_ = add_h2; }

   private:
    sigma_prior sigma_g_{4, 0};
    sigma_prior sigma_r_{2, 0};
    sigma_prior sigma_e_{-2, 0};
    dvec pi_;
    double add_h2_{0.5};
};
}  // namespace gelex
