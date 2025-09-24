#pragma once
#include <cstdint>
#include <optional>
#include <string_view>
#include <unordered_map>

#include <fmt/base.h>

namespace gelex
{

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
};

inline std::optional<BayesAlphabet> get_bayesalphabet(std::string_view sv)
{
    static const std::unordered_map<std::string_view, BayesAlphabet>
        stringToEnumMap = {
            {"A", BayesAlphabet::A},
            {"RR", BayesAlphabet::RR},
            {"B", BayesAlphabet::B},
            {"Bpi", BayesAlphabet::Bpi},
            {"C", BayesAlphabet::C},
            {"Cpi", BayesAlphabet::Cpi},
            {"R", BayesAlphabet::R},
        };

    auto it = stringToEnumMap.find(sv);

    if (it != stringToEnumMap.end())
    {
        return it->second;
    }
    return std::nullopt;
}

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
