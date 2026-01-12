#ifndef GELEX_TYPES_COVARIATES_H_
#define GELEX_TYPES_COVARIATES_H_

#include <string>
#include <vector>

#include <Eigen/Core>

namespace gelex
{

struct QuantitativeCovariate
{
    std::vector<std::string> names;
    Eigen::MatrixXd X;
};

struct DiscreteCovariate
{
    std::vector<std::string> names;
    std::vector<std::vector<std::string>> levels;
    std::vector<std::string> reference_levels;
    Eigen::MatrixXd X;
};

}  // namespace gelex

#endif  // GELEX_TYPES_COVARIATES_H_
