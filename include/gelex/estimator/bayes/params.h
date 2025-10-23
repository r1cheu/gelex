#pragma once
#include <cstddef>
#include <stdexcept>

#include <Eigen/Core>

namespace gelex
{
struct MCMCParams
{
    MCMCParams(
        Eigen::Index n_iters,
        Eigen::Index n_burnin,
        Eigen::Index n_thin,
        Eigen::Index n_chains)
        : n_iters{n_iters},
          n_burnin{n_burnin},
          n_thin{n_thin},
          n_chains{n_chains},
          n_records{(n_iters - n_burnin) / n_thin}
    {
        if (n_burnin >= n_iters)
        {
            throw std::invalid_argument(
                "n_burnin must be smaller than n_iters");
        }
    }
    Eigen::Index n_iters;
    Eigen::Index n_burnin;
    Eigen::Index n_thin;
    Eigen::Index n_chains;
    Eigen::Index n_records;
};
}  // namespace gelex
