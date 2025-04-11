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

double computeMeanDotProduct(const dmat& mat)
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

/**/
/*HyperParams::HyperParams(const BayesLMM& model)*/
/*{*/
/*    s2_ga_ = computeS2Ga(model);  // NOLINT*/
/*};*/
/**/
/*double HyperParams::computeS2Ga(const BayesLMM& model) const*/
/*{*/
/*    return ((nu_ga_ - 2) / nu_ga_) * model.var_y() / (1 - model.pi().at(0))*/
/*           * additive_h2_;*/
/*};*/
/**/
}  // namespace gelex
