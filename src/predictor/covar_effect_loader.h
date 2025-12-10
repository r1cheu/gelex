#ifndef GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H
#define GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H

#include <filesystem>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>

namespace gelex::detail
{

struct IndividualData
{
    std::unordered_map<std::string, double> continuous_values;
    std::unordered_map<std::string, std::string> categorical_values;
};

class CovarEffectLoader
{
   public:
    explicit CovarEffectLoader(const std::filesystem::path& param_file_path);

    double intercept() const { return intercept_; }
    const std::map<std::string, double>& continuous_coeffs() const
    {
        return continuous_coeffs_;
    }
    const std::map<std::string, std::map<std::string, double>>&
    categorical_coeffs() const
    {
        return categorical_coeffs_;
    }

    // Throws InvalidInputException if data contains variables or categories not
    // present in the parameter file
    double predict(const IndividualData& data) const;

   private:
    static auto parse_param_file(const std::filesystem::path& file_path)
        -> std::tuple<
            double,
            std::map<std::string, double>,
            std::map<std::string, std::map<std::string, double>>>;

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

#endif  // GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H
