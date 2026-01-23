#ifndef GELEX_ESTIMATOR_FREQ_REML_H_
#define GELEX_ESTIMATOR_FREQ_REML_H_

#include <memory>

#include "gelex/data/data_pipe.h"

namespace gelex
{
struct AssocInput;
class SampleManager;

auto load_data_for_reml(const DataPipe::Config& config) -> DataPipe;

auto reml(
    const DataPipe::Config& config,
    size_t max_iter = 100,
    double tol = 1e-8,
    bool em_init = true,
    bool verbose = true) -> std::
    tuple<std::shared_ptr<SampleManager>, Eigen::MatrixXd, Eigen::VectorXd>;

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_REML_H_
