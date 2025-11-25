#pragma once

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

namespace gelex
{

class PredictDataPipe;
struct PredictData
{
    explicit PredictData(PredictDataPipe& data_pipe);
    std::vector<std::string> individual_ids;

    Eigen::MatrixXd genotypes;

    Eigen::MatrixXd quantitative_covariates;
    std::vector<std::string> quantitative_covariate_names;

    std::map<std::string, std::vector<std::string>> categorical_covariates;

    void clear()
    {
        individual_ids.clear();
        genotypes.resize(0, 0);
        quantitative_covariates.resize(0, 0);
        quantitative_covariate_names.clear();
        categorical_covariates.clear();
    }
};

}  // namespace gelex
