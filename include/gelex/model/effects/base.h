#pragma once
#include <cstdint>

#include <fmt/base.h>
#include <armadillo>

namespace gelex
{
using dmat = arma::dmat;
using dvec = arma::dvec;

enum class effect_type : uint8_t
{
    random,
    genetic,
    gxe,
    residual,
};

enum class BayesAlphabet : uint8_t
{
    A,
    RR,
    B,
    Bpi,
    C,
    Cpi,
    R,
    None,
    Count
};

}  // namespace gelex

namespace fmt
{
template <>
struct formatter<gelex::BayesAlphabet> : formatter<string_view>
{
    auto format(gelex::BayesAlphabet, format_context& ctx) const
        -> format_context::iterator;
};

}  // namespace fmt
