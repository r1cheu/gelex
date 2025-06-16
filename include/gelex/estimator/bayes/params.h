#pragma once
#include <cstddef>
#include <stdexcept>

namespace gelex
{
struct MCMCParams
{
    MCMCParams(size_t n_iters, size_t n_burnin, size_t n_thin, size_t n_chains)
        : n_iters{n_iters},
          n_burnin{n_burnin},
          n_thin{n_thin},
          n_chains{n_chains}
    {
        if (n_burnin >= n_iters)
        {
            throw std::invalid_argument(
                "n_burnin must be smaller than n_iters");
        }
    }
    size_t n_iters;
    size_t n_burnin;
    size_t n_thin;
    size_t n_chains;
};
}  // namespace gelex
