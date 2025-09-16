#pragma once

#include <random>
#include <string>
#include <unordered_map>

#include <Eigen/Core>

namespace gelex
{

class PhenotypeSimulator
{
   public:
    PhenotypeSimulator();

    void simulate_qt_from_bed(
        const std::string& bfile,
        const std::string& causal_variants_list,
        double heritability,
        int seed = -1);

   private:
    std::mt19937_64 rng_;
    void initialize_rng(int seed);
    std::unordered_map<std::string, double> load_causal_variants(
        const std::string& causal_variants_list);
};

}  // namespace gelex
