#pragma once
#include <armadillo>
#include <random>

namespace gelex
{

arma::dvec dirichlet(const arma::uvec& alphas, std::mt19937_64& rng);

struct ScaledInvChiSqParams
{
    double nu{};
    double s2{};
};

inline double
sample_scale_inv_chi_squared(std::mt19937_64& rng, double nu, double s2)
{
    std::chi_squared_distribution<double> chisq{nu};
    return (nu * s2) / chisq(rng);
}

class ScaledInvChiSq
{
   public:
    explicit ScaledInvChiSq(const ScaledInvChiSqParams& prior_params);
    ScaledInvChiSq(double initial_nu, double initial_s2);

    void update(double sum_of_squared_errors, size_t num_observations);

    void update(double single_observation_squared_error);

    double sample(std::mt19937_64& rng) const;
    const ScaledInvChiSqParams& get_params() const;

   private:
    ScaledInvChiSqParams params_;
};
}  // namespace gelex
