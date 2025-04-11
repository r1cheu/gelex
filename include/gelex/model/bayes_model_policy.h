#pragma once
#include <armadillo>
#include <cstddef>

namespace gelex
{
using arma::dmat;
using arma::dvec;

template <typename T>
concept GeneticPolicy = requires {
    { T::init_pi() } -> std::same_as<dvec>;
} && requires {
    {
        T::init_sigma(std::declval<size_t>())
    } -> std::same_as<typename T::sigma_t>;
};

namespace detail
{
inline dvec default_sigma_vector(size_t n)
{
    return {n, arma::fill::zeros};
}
inline double default_sigma_scalar(size_t n)
{
    return 0;
}
inline dvec default_pi_2()
{
    return {0.95, 0.05};
}
inline dvec default_pi_4()
{
    return {0.95, 0.02, 0.02, 0.01};
}
}  // namespace detail

struct VectorSigmaPolicy
{
    using sigma_t = dvec;
    static dvec init_sigma(size_t n) { return detail::default_sigma_vector(n); }
};

struct ScalarSigmaPolicy
{
    using sigma_t = double;
    static double init_sigma(size_t n)
    {
        return detail::default_sigma_scalar(n);
    }
};

struct BayesAPolicy : VectorSigmaPolicy
{
    static dvec init_pi() { return {0.0, 1.0}; }
    static constexpr std::string name = "BayesA";
    static constexpr bool has_pi = false;
    static constexpr bool fixed_pi = false;
};

struct BayesBPolicy : VectorSigmaPolicy
{
    static dvec init_pi() { return detail::default_pi_2(); }
    static constexpr std::string name = "BayesB";
    static constexpr bool has_pi = true;
    static constexpr bool fixed_pi = false;
};

struct BayesRPolicy : VectorSigmaPolicy
{
    static dvec init_pi() { return detail::default_pi_4(); }
    static constexpr std::string name = "BayesR";
    static constexpr bool has_pi = true;
    static constexpr bool fixed_pi = true;
};

struct BayesBpiPolicy : VectorSigmaPolicy
{
    static dvec init_pi() { return detail::default_pi_2(); }
    static constexpr std::string name = "BayesBpi";
    static constexpr bool has_pi = true;
    static constexpr bool fixed_pi = true;
};

struct BayesCPolicy : ScalarSigmaPolicy
{
    static dvec init_pi() { return detail::default_pi_2(); }
    static constexpr std::string name = "BayesC";
    static constexpr bool has_pi = true;
    static constexpr bool fixed_pi = false;
};

struct BayesCpiPolicy : ScalarSigmaPolicy
{
    static dvec init_pi() { return detail::default_pi_2(); }
    static constexpr std::string name = "BayesCpi";
    static constexpr bool has_pi = true;
    static constexpr bool fixed_pi = true;
};

struct BayesRRPolicy : ScalarSigmaPolicy
{
    static dvec init_pi() { return {0.0, 1.0}; }
    static constexpr std::string name = "BayesRR";
    static constexpr bool has_pi = false;
    static constexpr bool fixed_pi = false;
};

static_assert(GeneticPolicy<BayesAPolicy>);
static_assert(GeneticPolicy<BayesBPolicy>);
static_assert(GeneticPolicy<BayesBpiPolicy>);
static_assert(GeneticPolicy<BayesCPolicy>);
static_assert(GeneticPolicy<BayesCpiPolicy>);
static_assert(GeneticPolicy<BayesRPolicy>);
static_assert(GeneticPolicy<BayesRRPolicy>);

}  // namespace gelex
