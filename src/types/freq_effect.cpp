#include "../src/types/freq_effect.h"

namespace gelex::freq
{

FixedEffect FixedEffect::build(
    QCovariateEffect&& qcovariate,
    DCovariateEffect&& dcovariate)
{
    QCovariateEffect qcov = std::move(qcovariate);
    DCovariateEffect dcov = std::move(dcovariate);
    FixedEffect fe;

    const auto n_samples = qcov.X.rows();
    const auto n_cols = 1 + qcov.X.cols() + dcov.X.cols();

    fe.names.reserve(n_cols);
    fe.names.emplace_back("Intercept");
    fe.names.insert(fe.names.end(), qcov.names.begin(), qcov.names.end());
    fe.names.insert(fe.names.end(), dcov.names.begin(), dcov.names.end());

    // QCovariateEffect does not have levels or reference levels
    fe.levels.assign(
        qcov.names.size() + 1, std::nullopt);  // add one for Intercept
    fe.levels.insert(fe.levels.end(), dcov.levels.begin(), dcov.levels.end());

    fe.reference_levels.assign(
        qcov.names.size() + 1, std::nullopt);  // add one for Intercept
    fe.reference_levels.insert(
        fe.reference_levels.end(),
        dcov.reference_levels.begin(),
        dcov.reference_levels.end());

    fe.X = Eigen::MatrixXd::Zero(n_samples, n_cols);
    fe.X.col(0).setOnes();
    fe.X.middleCols(1, qcov.X.cols()) = qcov.X;
    fe.X.middleCols(1 + qcov.X.cols(), dcov.X.cols()) = dcov.X;

    return fe;
}

}  // namespace gelex::freq
