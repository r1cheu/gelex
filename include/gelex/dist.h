#pragma once
#include <random>

#include <armadillo>

namespace gelex
{

using arma::dvec;
using arma::uvec;

dvec dirichlet(const uvec& alphas, std::mt19937_64& rng)
{
    dvec pi(alphas.size(), arma::fill::zeros);
    double sum = 0.0;
    for (size_t i = 0; i < alphas.size(); ++i)
    {
        std::gamma_distribution<double> gamma_dist(alphas[i], 1.0);
        pi.at(i) = gamma_dist(rng);
    }
    return pi / arma::sum(pi);
}

class Uniform
{
   public:
    explicit Uniform(std::mt19937_64& random_generator)
        : random_generator_(&random_generator), uniform_{0.0, 1.0}
    {
    }
    double operator()() noexcept { return uniform_(*random_generator_); }

   private:
    std::mt19937_64* random_generator_ = nullptr;
    std::uniform_real_distribution<double> uniform_;
};

class Normal
{
   public:
    explicit Normal(std::mt19937_64& random_generator)
        : random_generator_(&random_generator), normal_{0.0, 1.0}
    {
    }

    double operator()(double mu, double sigma) noexcept
    {
        return (normal_(*random_generator_) * sigma) + mu;
    }

   private:
    std::mt19937_64* random_generator_ = nullptr;
    std::normal_distribution<double> normal_;
};

class ScaleInvChiSq
{
   public:
    ScaleInvChiSq(
        std::mt19937_64& random_generator,
        double prior_nu,
        double nu,
        double s2)
        : random_generator_(&random_generator),
          prior_nu_{prior_nu},
          chi_squared_{prior_nu + nu},
          s2_adj_{s2 * prior_nu}
    {
    }

    double update(double nu, double ssq)
    {
        chi_squared_ = std::chi_squared_distribution<double>(prior_nu_ + nu);
        return (ssq + s2_adj_) / chi_squared_(*random_generator_);
    }

    double operator()(double ssq) noexcept
    {
        return (ssq + s2_adj_) / chi_squared_(*random_generator_);
    }

   private:
    std::mt19937_64* random_generator_ = nullptr;
    double prior_nu_{};
    double s2_adj_{};
    std::chi_squared_distribution<double> chi_squared_;
};
}  // namespace gelex
