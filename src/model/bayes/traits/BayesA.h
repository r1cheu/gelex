#pragma once
#include "base_trait.h"

#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"

namespace gelex
{
class BayesATrait : public GeneticTrait
{
   public:
    std::string name() const override { return "BayesA"; }

    Eigen::VectorXd default_marker_variance(Eigen::Index n_snp) const override;
    Eigen::VectorXd default_pi() const override;
    bool estimate_pi() const override;

    std::vector<std::string> prior_info(
        double nu,
        double s2,
        const Eigen::Ref<const Eigen::VectorXd>& pi) const override;

    void operator()(
        const bayes::AdditiveEffect& effect,
        bayes::AdditiveStatus& state,
        Eigen::Ref<Eigen::VectorXd> y_adj,
        double sigma_e,
        std::mt19937_64& rng) const override;
};
}  // namespace gelex
