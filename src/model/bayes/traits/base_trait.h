#pragma once

#include <Eigen/Core>
#include <vector>

#include "../src/model/bayes/bayes_effects.h"

namespace gelex
{

class GeneticTrait
{
   public:
    virtual ~GeneticTrait() = default;

    virtual Eigen::VectorXd default_marker_variance(Eigen::Index n_snp) const
        = 0;
    virtual Eigen::VectorXd default_pi() const = 0;
    virtual bool estimate_pi() const = 0;

    virtual std::vector<std::string> prior_info(
        double nu,
        double s2,
        const Eigen::Ref<const Eigen::VectorXd>& pi) const
        = 0;

    virtual void operator()(
        const bayes::AdditiveEffect& effect,
        bayes::AdditiveStatus& state,
        Eigen::Ref<Eigen::VectorXd> y_adj,
        double sigma_e,
        std::mt19937_64& rng) const
        = 0;

    virtual std::string name() const = 0;
};

}  // namespace gelex
