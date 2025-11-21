#pragma once

#include <Eigen/Dense>
#include <map>
#include <string>
#include <vector>

namespace gelex
{

struct PredictData
{
    std::vector<std::string> individual_ids;

    Eigen::MatrixXd genotypes;

    Eigen::MatrixXd quantitative_covariates;
    std::vector<std::string> quantitative_covariate_names;

    // key: covariate name, value: categorical values for each individual
    std::map<std::string, std::vector<std::string>> categorical_covariates;

    std::string bed_file_path;
    std::string bim_file_path;
    std::string fam_file_path;

    void clear()
    {
        individual_ids.clear();
        genotypes.resize(0, 0);
        quantitative_covariates.resize(0, 0);
        quantitative_covariate_names.clear();
        categorical_covariates.clear();
        bed_file_path.clear();
        bim_file_path.clear();
        fam_file_path.clear();
    }
};

struct PredictResult
{
    std::vector<std::string> individual_ids;

    Eigen::VectorXd genetic_values;
    Eigen::VectorXd covariate_effects;
    Eigen::VectorXd total_predictions;

    void clear()
    {
        individual_ids.clear();
        genetic_values.resize(0);
        covariate_effects.resize(0);
        total_predictions.resize(0);
    }
};

}  // namespace gelex
