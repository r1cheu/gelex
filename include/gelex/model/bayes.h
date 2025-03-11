#include <armadillo>
#include <random>

namespace gelex
{
using arma::dmat;
using arma::dvec;
static constexpr uint64_t DEFAULT_N_ITER = 20000;
static constexpr uint64_t DEFAULT_N_BURN = 2000;
static constexpr uint64_t SEED = 42;

class BayesRR
{
   public:
    BayesRR(
        dvec&& phenotype,
        dmat&& design_matrix_beta,
        dmat&& genotype_matrix,
        uint64_t seed = SEED,
        uint64_t n_iter = DEFAULT_N_ITER,
        uint64_t n_burn = DEFAULT_N_BURN);
    void set_hyperparameters(double nu_e, double nu_a);

   private:
    // --- functions ---
    void init_params();
    double SampleElement(
        const dvec& col_vector,
        double col_norm2,
        double current_effect,
        double sigma_e,
        double regularization_ratio);
    void SamplingBeta(uint64_t index);
    void SamplingSnpEffect(uint64_t index);

    // --- data ---
    dvec phenotype_;
    dmat genotype_matrix_;
    dmat design_matrix_beta_;
    std::mt19937 rng_;
    std::normal_distribution<double> normal_{0.0, 1.0};

    double rss_{};

    dvec b_cols_norm2_;  // design_matrix_beta_.col(index).T *
                         // design_matrix_beta_.col(index)
    dvec g_cols_norm2_;  // genotype_matrix_.col(index).T *
                         // genotype_matrix_.col(index)

    uint64_t n_iter_;
    uint64_t n_burn_;

    // --- parameters ---
    double mu_{};
    dvec beta_;
    dvec snp_effect_;

    // --- mcmc sample sum --
    double mu_sum_{};
    dvec beta_sum_;
    dvec snp_effect_sum_;

    double sigma_e_{};
    double sigma_a_{};

    // --- hyper-paramters ---
    double nu_e_{-2};
    double s_e0_{};
    double nu_a_{4};
    double s_a0_{};
};

}  // namespace gelex
