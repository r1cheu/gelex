/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef GELEX_PREDICT_COVAR_EFFECT_LOADER_H
#define GELEX_PREDICT_COVAR_EFFECT_LOADER_H

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

#endif  // GELEX_PREDICT_COVAR_EFFECT_LOADER_H
