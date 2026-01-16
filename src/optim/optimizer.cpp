#include "gelex/optim/optimizer.h"

namespace gelex
{

auto collect_variance_components(const FreqState& state) -> Eigen::VectorXd
{
    auto n_random = static_cast<Eigen::Index>(state.random().size());
    auto n_genetic = static_cast<Eigen::Index>(state.genetic().size());
    Eigen::Index n_total
        = 1 + n_random + n_genetic;  // residual + random + genetic

    Eigen::VectorXd sigma(n_total);

    // residual variance first
    sigma(0) = state.residual().variance;

    // random effects
    Eigen::Index idx = 1;
    for (const auto& r : state.random())
    {
        sigma(idx++) = r.variance;
    }

    // genetic effects
    for (const auto& g : state.genetic())
    {
        sigma(idx++) = g.variance;
    }

    return sigma;
}

auto distribute_variance_components(
    FreqState& state,
    const Eigen::Ref<const Eigen::VectorXd>& sigma) -> void
{
    // residual variance first
    state.residual().variance = sigma(0);

    // random effects
    Eigen::Index idx = 1;
    for (auto& r : state.random())
    {
        r.variance = sigma(idx++);
    }

    // genetic effects
    for (auto& g : state.genetic())
    {
        g.variance = sigma(idx++);
    }
}

}  // namespace gelex
