#include "gelex/model/effects/base.h"
#include "fmt/format.h"

namespace fmt
{

auto formatter<gelex::BayesAlphabet>::format(
    gelex::BayesAlphabet t,
    format_context& ctx) const -> format_context::iterator
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
            name = "None";
            break;
    }
    return formatter<string_view>::format(name, ctx);
}
}  // namespace fmt
