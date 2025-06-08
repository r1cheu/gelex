#pragma once
#include <cstddef>

namespace gelex
{
struct MCMCParams
{
    size_t iter;
    size_t n_burnin;
    size_t n_thin;
    size_t n_chains;
};
}  // namespace gelex
