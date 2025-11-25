#pragma once

#include <vector>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"

namespace gelex
{
struct MatchInfo;

class PredictBedPipe
{
   public:
    PredictBedPipe(
        const std::filesystem::path& bed_path,
        const std::filesystem::path& snp_effect_path,
        std::shared_ptr<SampleManager> sample_manager);

    auto load() const -> std::vector<Eigen::MatrixXd>;

   private:
    BedPipe bed_pipe_;
    bool has_dom_ = false;
};

}  // namespace gelex
