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
    auto format(gelex::BayesAlphabet t, format_context& ctx) const
        -> format_context::iterator
    {
        string_view name = "unknown";
        switch (t)
        {
            case gelex::BayesAlphabet::A:
                name = "BayesA";
                break;
            case gelex::BayesAlphabet::RR:
                name = "BayesRR";
                break;
            case gelex::BayesAlphabet::B:
                name = "BayesB";
                break;
            case gelex::BayesAlphabet::Bpi:
                name = "BayesBpi";
                break;
            case gelex::BayesAlphabet::C:
                name = "BayesC";
                break;
            case gelex::BayesAlphabet::Cpi:
                name = "BayesCpi";
                break;
            case gelex::BayesAlphabet::R:
                name = "BayesR";
                break;
            default:
                break;
        }
        return formatter<string_view>::format(name, ctx);
    }
};

}  // namespace fmt
