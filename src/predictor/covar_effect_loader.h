#ifndef GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H
#define GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H

#include <filesystem>
#include <limits>
#include <map>
#include <string>

namespace gelex::detail
{

struct CovarEffects
{
    double intercept = std::numeric_limits<double>::quiet_NaN();
    std::map<std::string, double> continuous_coeffs;
    std::map<std::string, std::map<std::string, double>> categorical_coeffs;
};

class CovarEffectLoader
{
   public:
    explicit CovarEffectLoader(const std::filesystem::path& param_file_path);
    const CovarEffects& effects() const { return effects_; }
    CovarEffects&& take_effects() && { return std::move(effects_); }

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

    CovarEffects effects_;
};

}  // namespace gelex::detail

#endif  // GELEX_PREDICTOR_COVAR_EFFECT_LOADER_H
