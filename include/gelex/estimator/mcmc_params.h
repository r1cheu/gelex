#pragma once
#include <cstddef>

namespace gelex
{
struct MCMCParams
{
    size_t iter;
    size_t n_burnin;
    size_t n_thin;
    size_t seed;

    MCMCParams(size_t iter, size_t n_burnin, size_t n_thin, size_t seed)
        : iter{iter}, n_burnin{n_burnin}, n_thin{n_thin}, seed{seed}
    {
    }
};
}  // namespace gelex
