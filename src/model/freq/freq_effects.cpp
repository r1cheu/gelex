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
}  // namespace freq

}  // namespace gelex
