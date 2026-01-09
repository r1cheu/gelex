#pragma once
#include <string>
#include <vector>

#include <armadillo>

#include "../effects_manager.h"

namespace gelex
{

namespace freq
{
struct FixedEffect
{
    FixedEffect(
        std::vector<std::string>&& names,
        std::vector<std::string>&& levels,
        arma::dmat&& design_matrix);

    size_t size() const { return coeff.n_elem; }

    std::vector<std::string> names;
    std::vector<std::string> levels;

    arma::dmat design_matrix;
    arma::dvec coeff;
};

struct RandomEffect
{
    RandomEffect(std::string&& name, arma::sp_dmat&& design_matrix);

    std::string name;

    arma::sp_dmat design_matrix;
    arma::sp_dmat covariance_matrix;

    arma::dvec coeff;
    double sigma{};
};

struct GeneticEffect
{
    GeneticEffect(
        std::string&& name,
        arma::sp_dmat&& design_matrix,
        const arma::dmat& genetic_relationship_matrix);

    std::string name;

    arma::sp_dmat design_matrix;
    arma::dmat genetic_relationship_matrix;
    arma::dmat covariance_matrix;

    arma::dvec coeff;
    double sigma{};
};

using RandomEffects = detail::Effects<RandomEffect>;
using FixedEffects = detail::Effects<FixedEffect>;
using GeneticEffects = detail::Effects<GeneticEffect>;

}  // namespace freq
}  // namespace gelex
