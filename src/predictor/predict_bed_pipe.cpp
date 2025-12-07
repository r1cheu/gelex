#include "../src/predictor/predict_bed_pipe.h"

#include "../src/predictor/snp_matcher.h"
#include "gelex/data/bed_pipe.h"

#include <limits>

namespace gelex
{

PredictBedPipe::PredictBedPipe(
    const std::filesystem::path& bed_path,
    const std::filesystem::path& snp_effect_path,
    std::shared_ptr<SampleManager> sample_manager)
    : bed_pipe_(bed_path, std::move(sample_manager)),
      has_dom_(detail::has_dom_effect_column(snp_effect_path))
{
    SnpMatcher matcher(snp_effect_path);
    match_plan_ = matcher.match(bed_path);
    snp_effects_ = std::move(matcher).take_snp_effects();
}

auto PredictBedPipe::load() const -> std::vector<Eigen::MatrixXd>
{
    const Eigen::MatrixXd full_matrix = bed_pipe_.load();
    const Eigen::Index num_samples = bed_pipe_.num_samples();
    const Eigen::Index num_effect_snps = snp_effects_.size();

    Eigen::MatrixXd add_matrix(num_samples, num_effect_snps);
    add_matrix.setConstant(std::numeric_limits<double>::quiet_NaN());

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

        add_matrix.col(info.target_col) = col;
    }

    std::vector<Eigen::MatrixXd> genotypes;
    genotypes.push_back(std::move(add_matrix));

    if (has_dom_)
    {
        Eigen::MatrixXd dom_matrix(num_samples, num_effect_snps);
        dom_matrix.setZero();
        genotypes.push_back(std::move(dom_matrix));
    }

    return genotypes;
}
}  // namespace gelex
