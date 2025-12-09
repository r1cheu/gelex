#include "../src/predictor/predict_bed_pipe.h"

#include <Eigen/Core>

#include "../src/predictor/snp_matcher.h"
#include "gelex/data/bed_pipe.h"

namespace gelex
{

PredictBedPipe::PredictBedPipe(
    const std::filesystem::path& bed_path,
    const std::filesystem::path& snp_effect_path,
    std::shared_ptr<SampleManager> sample_manager)
    : bed_pipe_(bed_path, std::move(sample_manager))
{
    SnpMatcher matcher(snp_effect_path);
    match_plan_ = matcher.match(bed_path);
    snp_effects_ = std::move(matcher).take_snp_effects();
}

auto PredictBedPipe::load() const -> Eigen::MatrixXd
{
    const Eigen::MatrixXd full_matrix = bed_pipe_.load();
    const Eigen::Index num_samples = bed_pipe_.num_samples();
    const auto num_effect_snps = static_cast<Eigen::Index>(snp_effects_.size());

    Eigen::MatrixXd genotype(num_samples, num_effect_snps);
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
