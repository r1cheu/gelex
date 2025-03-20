#include "gelex/model/bayes_prior.h"

#include <armadillo>
#include <cmath>

namespace gelex
{
using arma::dmat;
using arma::dvec;

double SumColVar(const dmat& mat)
{
    dvec out(mat.n_cols);
#pragma omp parallel for default(none) shared(mat, out)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        out.at(i) = arma::var(mat.col(i));
    }
    return arma::sum(out);
}

double ComputeMeanDotProduct(const dmat& mat)
{
    double sum = 0.0;
#pragma omp parallel for default(none) shared(mat) reduction(+ : sum) \
    schedule(static)
    for (size_t i = 0; i < mat.n_cols; ++i)
    {
        double col_mean = arma::mean(mat.col(i));
        sum += 2 * col_mean * (1 - col_mean);
    }
    return sum;
}

void SimgaPriors::set_genomic_effect_s2(const dvec& y, const dmat& genotype_mat)
{
    if (std::fabs(sigma_g_.s2) < 1e-8)  // check if user has set s2
    {
        return;
    }
    const double scale_term{
        (1 - pi_.at(0)) * ComputeMeanDotProduct(genotype_mat)};
    const double y_var{arma::var(y)};

    sigma_g_.s2 = ((sigma_g_.nu - 2) / sigma_g_.nu) * y_var / (1 - pi_.at(0))
                  / scale_term * add_h2_;
}

/**/
/*HyperParams::HyperParams(const BayesLMM& model)*/
/*{*/
/*    s2_ga_ = ComputeS2Ga(model);  // NOLINT*/
/*};*/
/**/
/*double HyperParams::ComputeS2Ga(const BayesLMM& model) const*/
/*{*/
/*    return ((nu_ga_ - 2) / nu_ga_) * model.var_y() / (1 - model.pi().at(0))*/
/*           * additive_h2_;*/
/*};*/
/**/
}  // namespace gelex
