#include "../src/predict/predict_bed_pipe.h"

#include <Eigen/Core>

#include "../src/predict/snp_matcher.h"
#include "gelex/data/bed_pipe.h"

namespace gelex
{

PredictBedPipe::PredictBedPipe(
    const std::filesystem::path& bed_path,
    const SnpEffects& snp_effects,
    std::shared_ptr<SampleManager> sample_manager)
    : bed_pipe_(bed_path, std::move(sample_manager)),
      num_snp_effects_(static_cast<Eigen::Index>(snp_effects.size()))
{
    SnpMatcher matcher(snp_effects);
    match_plan_ = matcher.match(bed_path);
}

auto PredictBedPipe::load() const -> Eigen::MatrixXd
{
    const Eigen::MatrixXd full_matrix = bed_pipe_.load();
    const Eigen::Index num_samples = bed_pipe_.num_samples();

    Eigen::MatrixXd genotype(num_samples, num_snp_effects_);
    genotype.setZero();

#pragma omp parallel for schedule(dynamic) default(none) \
    shared(full_matrix, genotype)
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(match_plan_.size());
         ++i)
    {
        const MatchInfo& info = match_plan_[i];

        if (info.type == MatchType::skip || info.target_col < 0)
        {
            continue;
        }

        Eigen::VectorXd col = full_matrix.col(i);

        if (info.type == MatchType::reverse)
        {
            col = 2.0 - col.array();
        }

        genotype.col(info.target_col) = col;
    }

    return genotype;
}
}  // namespace gelex
