#pragma once
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <armadillo>

namespace gelex
{
using arma::dcube;
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
using arma::uvec;
using arma::uword;

class GBLUP
{
   public:
    GBLUP(dvec phenotype, dmat design_mat_beta);

    uint64_t n_individuals() const { return n_individuals_; }
    uint64_t n_common_effects() const { return n_common_effects_; }
    uint64_t n_group_effects() const { return n_group_effects_; }
    uint64_t n_genetic_effects() const { return n_genetic_effects_; }
    uint64_t n_random_effects() const { return n_random_effects_; }

    const dvec& phenotype() const { return phenotype_; }
    const dmat& design_mat_beta() const { return design_mat_beta_; }
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }
    const std::vector<std::string>& sigma_names() const { return sigma_names_; }

    void add_group_effect(std::string name, sp_dmat design_mat_env);

    void add_genetic_effect(std::string name, dmat genetic_covar_mat);

    const std::vector<dmat>& genetic_cov_mats() const
    {
        return genetic_cov_mats_;
    }

    const std::vector<sp_dmat>& env_cov_mats() const { return env_cov_mats_; }

    void set_sigma(dvec sigma) { sigma_ = std::move(sigma); }
    void set_beta(dvec beta) { beta_ = std::move(beta); }

    void set_model();
    void reset();

   private:
    uint64_t n_individuals_{};
    uint64_t n_common_effects_{};
    uint64_t n_group_effects_{};
    uint64_t n_genetic_effects_{};
    uint64_t n_random_effects_{};

    dvec phenotype_;
    dmat design_mat_beta_;

    std::vector<sp_dmat> env_cov_mats_;
    std::vector<dmat> genetic_cov_mats_;

    std::vector<std::string> sigma_names_;

    dvec beta_;
    dvec sigma_;
};

class GBLUPParams
{
   public:
    GBLUPParams(
        dvec beta,
        dvec sigma,
        dvec proj_y,
        std::vector<std::string> dropped_individuals);
    const dvec& beta() const { return beta_; }
    const dvec& sigma() const { return sigma_; }
    const dvec& proj_y() const { return proj_y_; }
    const std::vector<std::string>& dropped_individuals() const
    {
        return dropped_individuals_;
    }

    void set_beta(dvec beta) { beta_ = std::move(beta); }
    void set_sigma(dvec sigma) { sigma_ = std::move(sigma); }
    void set_proj_y(dvec proj_y) { proj_y_ = std::move(proj_y); }
    void set_dropped_individuals(std::vector<std::string> dropped_individuals)
    {
        dropped_individuals_ = std::move(dropped_individuals);
    }

   private:
    dvec beta_;
    dvec sigma_;
    dmat X_;
    dvec proj_y_;
    std::vector<std::string> dropped_individuals_;
};

}  // namespace gelex
