#include "../src/predictor/predict_bed_pipe.h"

#include "../src/predictor/snp_matcher.h"
#include "gelex/data/bed_pipe.h"

namespace gelex
{

PredictBedPipe::PredictBedPipe(
    const std::filesystem::path& bed_path,
    const std::filesystem::path& snp_effect_path,
    std::shared_ptr<SampleManager> sample_manager)
    : bed_pipe_(bed_path, std::move(sample_manager)),
      has_dom_(SnpEffectLoader::has_dom_effects(snp_effect_path))
{
    SnpMatcher matcher(snp_effect_path);
    MatchPlan plan = matcher.match(bed_path);

    bed_pipe_.set_read_plan(std::move(plan));
}

auto PredictBedPipe::load() const -> std::vector<Eigen::MatrixXd>
{
    std::vector<Eigen::MatrixXd> genotypes;
    genotypes.push_back(bed_pipe_.load());
    if (has_dom_)
    {
        // TODO(rlchen): fill in dom genotype loading
    }
    return genotypes;
}
}  // namespace gelex
