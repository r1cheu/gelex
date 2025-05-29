#pragma once
#include <cstddef>

namespace gelex
{
struct MCMCParams
{
    size_t iter;
    size_t n_burnin;
    size_t n_thin;

    MCMCParams(size_t iter, size_t n_burnin, size_t n_thin)
        : iter{iter}, n_burnin{n_burnin}, n_thin{n_thin}
    {
    }
};
}  // namespace gelex
