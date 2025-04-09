#pragma once
#include <armadillo>

namespace gelex
{
using arma::dmat;
using arma::dvec;
using arma::uvec;

struct sigma_prior
{
    double nu;
    double s2;
};

class Priors
{
   public:
    Priors() = default;
    explicit Priors(const dvec& pi) : pi_{pi} {}

    sigma_prior& sigma_a() { return sigma_a_; }
    sigma_prior& sigma_r() { return sigma_r_; }
    sigma_prior& sigma_e() { return sigma_e_; }

    dvec pi() const { return pi_; }
    void set_pi(dvec pi) { pi_ = std::move(pi); }

   private:
    sigma_prior sigma_a_{-2, 0};
    sigma_prior sigma_r_{2, 0};
    sigma_prior sigma_e_{-2, 0};
    dvec pi_;
};
}  // namespace gelex
