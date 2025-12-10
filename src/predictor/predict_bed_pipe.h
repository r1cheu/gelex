#ifndef GELEX_PREDICTOR_PREDICT_BED_PIPE_H
#define GELEX_PREDICTOR_PREDICT_BED_PIPE_H

#include <Eigen/Core>

#include "../src/predictor/snp_matcher.h"
#include "gelex/data/bed_pipe.h"

namespace gelex
{
class PredictBedPipe
{
   public:
    PredictBedPipe(
        const std::filesystem::path& bed_path,
        const SnpEffects& snp_effects,
        std::shared_ptr<SampleManager> sample_manager);

    auto load() const -> Eigen::MatrixXd;

   private:
    BedPipe bed_pipe_;
    MatchPlan match_plan_;
    Eigen::Index num_snp_effects_ = 0;
};

}  // namespace gelex

#endif  // GELEX_PREDICTOR_PREDICT_BED_PIPE_H
