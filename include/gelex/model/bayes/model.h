#pragma once

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "../src/model/bayes/bayes_effects.h"
#include "../src/model/bayes/trait.h"
#include "gelex/model/effects.h"

namespace gelex
{
/**
 * @class BayesModel
 * @brief Bayesian model for genetic analysis with fixed, random, and genetic
 * effects.
 *
 * Supports adding and managing effects (fixed, random, genetic) and provides
 * accessors for model components (mu, residuals, etc.). Initialized with a
 * formula and phenotype.
 */
class BayesModel
{
   public:
    /**
     * @brief
     *
     * @param formula formula string, only for showing the model structure.
     * @param phenotype moveable arma dvec.
     */
    explicit BayesModel(Eigen::VectorXd&& phenotype);
    void add_fixed_effect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        Eigen::MatrixXd&& design_matrix);
    void add_random_effect(std::string&& name, Eigen::MatrixXd&& design_matrix);
    void add_genetic_effect(
        bool iid_only,
        bool dom,
        const std::string& bfile,
        const std::vector<std::string>& pid,
        const std::vector<std::string>& gid,
        BayesAlphabet type);

    const auto& fixed() const { return *fixed_; }
    const auto& random() const { return random_; }
    const auto& additive() const { return additive_; }
    const auto& dominant() const { return dominant_; }
    const auto& residual() const { return residual_; }

    auto& fixed() { return *fixed_; }
    auto& random() { return random_; }
    auto& additive() { return additive_; }
    auto& dominant() { return dominant_; }
    auto& residual() { return residual_; }

    void set_model_type(BayesAlphabet type)
    {
        model_trait_ = create_genetic_trait(type);
    }

    void set_sigma_prior_manual(
        const std::string& name,
        double nu,
        double s2,
        double init_sigma);

    void set_sigma_prior(const std::string& name, double nu, double prop);

    const std::unique_ptr<GeneticTrait>& trait() const { return model_trait_; };

    // Keep existing pi prior method
    void set_pi_prior(const std::string& name, const Eigen::VectorXd& pi);

    std::string prior_summary() const;

    const Eigen::VectorXd& phenotype() const { return phenotype_; }
    double phenotype_var() const { return phenotype_var_; }
    Eigen::Index n_individuals() const { return n_individuals_; }

   private:
    Eigen::Index n_individuals_{};

    Eigen::VectorXd phenotype_;
    double phenotype_var_{};

    std::unique_ptr<GeneticTrait> model_trait_;

    std::unique_ptr<bayes::FixedEffect> fixed_;
    bayes::RandomEffectManager random_;
    std::unique_ptr<bayes::AdditiveEffect> additive_;
    std::unique_ptr<bayes::DominantEffect> dominant_;
    bayes::Residual residual_;
};

struct BayesStatus
{
    explicit BayesStatus(const BayesModel& model)
        : fixed(bayes::FixedStatus(model.fixed().design_matrix.cols())),
          random(bayes::create_chain_states(model.random())),
          additive(
              model.additive()
                  ? std::make_unique<bayes::AdditiveStatus>(*model.additive())
                  : nullptr),
          dominant(
              model.dominant()
                  ? std::make_unique<bayes::DominantStatus>(*model.dominant())
                  : nullptr),
          residual(model.residual())
    {
    }

    void compute_heritability();

    bayes::FixedStatus fixed;
    std::vector<bayes::RandomStatus> random;
    std::unique_ptr<bayes::AdditiveStatus> additive;
    std::unique_ptr<bayes::DominantStatus> dominant;
    bayes::Residual residual;
};

}  // namespace gelex
