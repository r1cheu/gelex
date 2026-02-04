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

#ifndef GELEX_DATA_SIMULATE_H_
#define GELEX_DATA_SIMULATE_H_

#include <filesystem>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>

namespace gelex
{
class BedPipe;
class SampleManager;

struct EffectSizeClass
{
    double proportion{1.0};
    double variance{1.0};
};

struct CausalEffect
{
    double additive{0.0};
    double dominance{0.0};
    int add_class{0};
    int dom_class{0};
};

struct GeneticValues
{
    Eigen::VectorXd additive;
    Eigen::VectorXd dominance;
};

class PhenotypeSimulator
{
   public:
    struct Config
    {
        std::filesystem::path bed_path;
        double add_heritability{0.5};
        double dom_heritability{0.0};
        std::vector<EffectSizeClass> add_effect_classes{{1.0, 1.0}};
        std::vector<EffectSizeClass> dom_effect_classes{{1.0, 1.0}};
        double intercept{0.0};
        int seed{-1};
        std::filesystem::path output_path;
    };

    explicit PhenotypeSimulator(Config config);

    PhenotypeSimulator(const PhenotypeSimulator&) = delete;
    PhenotypeSimulator(PhenotypeSimulator&&) noexcept = default;
    PhenotypeSimulator& operator=(const PhenotypeSimulator&) = delete;
    PhenotypeSimulator& operator=(PhenotypeSimulator&&) noexcept = default;
    ~PhenotypeSimulator() = default;

    void simulate();

   private:
    void initialize_rng();

    auto select_causal_snps(const std::vector<std::string>& snp_ids)
        -> std::unordered_map<Eigen::Index, CausalEffect>;

    auto calculate_genetic_values(
        BedPipe& bed_pipe,
        const std::unordered_map<Eigen::Index, CausalEffect>& causal_effects)
        const -> GeneticValues;

    auto generate_phenotypes(
        const Eigen::VectorXd& additive_values,
        const Eigen::VectorXd& dominance_values) -> Eigen::VectorXd;

    void write_results(
        const Eigen::VectorXd& phenotypes,
        const std::shared_ptr<SampleManager>& sample_manager) const;

    void write_causal_effects(
        const std::vector<std::string>& snp_ids,
        const std::unordered_map<Eigen::Index, CausalEffect>& causal_effects)
        const;

    void write_params() const;

    Config config_;
    std::mt19937_64 rng_;
    double true_h2_{0.0};
    double true_d2_{0.0};

    static constexpr Eigen::Index SNP_CHUNK_SIZE = 10000;
};

}  // namespace gelex

#endif  // GELEX_DATA_SIMULATE_H_
