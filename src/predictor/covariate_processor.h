#pragma once

#include <expected>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>

#include "gelex/error.h"

namespace gelex::detail
{

struct IndividualData
{
    std::unordered_map<std::string, double> continuous_values;
    std::unordered_map<std::string, std::string> categorical_values;
};

class CovariateProcessor
{
   public:
    static auto create(const std::filesystem::path& param_file_path)
        -> std::expected<CovariateProcessor, Error>;

    double predict(const IndividualData& data) const;

   private:
    explicit CovariateProcessor(
        double intercept,
        std::map<std::string, double> continuous_coeffs,
        std::map<std::string, std::map<std::string, double>>
            categorical_coeffs);

    static auto parse_param_file(const std::filesystem::path& file_path)
        -> std::expected<
            std::tuple<
                double,
                std::map<std::string, double>,
                std::map<std::string, std::map<std::string, double>>>,
            Error>;

    static void parse_flat_name(
        const std::string& flat_name,
        double coefficient,
        double& intercept,
        std::map<std::string, double>& continuous_coeffs,
        std::map<std::string, std::map<std::string, double>>&
            categorical_coeffs);

    double intercept_ = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> continuous_coeffs_;
    std::map<std::string, std::map<std::string, double>> categorical_coeffs_;
};

}  // namespace gelex::detail
