#ifndef GELEX_PREDICT_GENOTYPE_ALIGNER_H
#define GELEX_PREDICT_GENOTYPE_ALIGNER_H

#include <Eigen/Core>

#include "../src/predict/snp_matcher.h"

namespace gelex
{
class GenotypeAligner
{
   public:
    GenotypeAligner(
        const std::filesystem::path& bed_path,
        const SnpEffects& snp_effects);

    auto load(Eigen::MatrixXd&& raw_genotype) const -> Eigen::MatrixXd;

   private:
    MatchPlan match_plan_;
    Eigen::Index num_snp_effects_ = 0;
};

}  // namespace gelex

#endif  // GELEX_PREDICT_GENOTYPE_ALIGNER_H
