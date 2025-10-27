#include <string>
#include <vector>

#include <armadillo>
#include "../src/utils/utils.h"

#include "freq_effects.h"

namespace gelex
{
namespace freq
{

FixedEffect::FixedEffect(
    std::vector<std::string>&& names,
    std::vector<std::string>&& levels,
    arma::dmat&& design_matrix)
    : names(std::move(names)),
      levels(std::move(levels)),
      design_matrix(std::move(design_matrix))
{
    coeff = arma::zeros(this->design_matrix.n_cols);
}

RandomEffect::RandomEffect(std::string&& name, arma::sp_dmat&& design_matrix)
    : name(std::move(name)), design_matrix(std::move(design_matrix))
{
    if (check_eye(this->design_matrix))
    {
        covariance_matrix = this->design_matrix;
    }
    else
    {
        covariance_matrix = this->design_matrix * this->design_matrix.t();
    }
    coeff = arma::zeros(this->design_matrix.n_cols);
}
GeneticEffect::GeneticEffect(
    std::string&& name,
    arma::sp_dmat&& design_matrix,
    const arma::dmat& genetic_relationship_matrix)

    : name(std::move(name)),
      design_matrix(std::move(design_matrix)),
      genetic_relationship_matrix(genetic_relationship_matrix)
{
    if (check_eye(this->design_matrix))
    {
        covariance_matrix = this->genetic_relationship_matrix;
    }
    else
    {
        covariance_matrix = this->design_matrix
                            * this->genetic_relationship_matrix
                            * this->design_matrix.t();
    }
    coeff = arma::zeros(this->design_matrix.n_cols);
}

GxEEffect::GxEEffect(
    std::string&& name,
    arma::sp_dmat&& genetic_design_matrix,
    const arma::dmat& genetic_relationship_matrix,
    arma::sp_dmat&& env_design_matrix)
    : name(std::move(name)),
      genetic_design_matrix(std::move(genetic_design_matrix)),
      genetic_relationship_matrix(genetic_relationship_matrix),
      env_design_matrix(std::move(env_design_matrix))
{
    arma::sp_dmat cov_env
        = this->env_design_matrix * this->env_design_matrix.t();

    if (check_eye(this->genetic_design_matrix))
    {
        covariance_matrix = this->genetic_relationship_matrix % cov_env;
    }
    else
    {
        covariance_matrix = this->genetic_design_matrix
                            * this->genetic_relationship_matrix
                            * this->genetic_design_matrix.t();
        covariance_matrix %= cov_env;
    }
    coeff = arma::zeros(this->genetic_design_matrix.n_cols);
}
}  // namespace freq

}  // namespace gelex
