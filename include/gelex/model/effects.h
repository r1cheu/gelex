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
    Ad,
    RRd,
    Bd,
    Bdpi,
    Cd,
    Cdpi,
    Rd,
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
            {"Ad", BayesAlphabet::Ad},
            {"RRd", BayesAlphabet::RRd},
            {"Bd", BayesAlphabet::Bd},
            {"Bdpi", BayesAlphabet::Bdpi},
            {"Cd", BayesAlphabet::Cd},
            {"Cdpi", BayesAlphabet::Cdpi},
            {"Rd", BayesAlphabet::Rd},
        };

    auto it = stringToEnumMap.find(sv);

    if (it != stringToEnumMap.end())
    {
        return it->second;
    }
    return std::nullopt;
}

}  // namespace gelex
