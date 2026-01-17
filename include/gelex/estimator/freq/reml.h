#ifndef GELEX_ESTIMATOR_FREQ_REML_H_
#define GELEX_ESTIMATOR_FREQ_REML_H_

#include <memory>
#include <utility>

#include "gelex/data/data_pipe.h"
#include "gelex/types/assoc_input.h"

namespace gelex
{
struct AssocInput;
class SampleManager;

auto reml(
    const DataPipe::Config& config,
    size_t max_iter = 100,
    double tol = 1e-8,
    bool em_init = true,
    bool verbose = true)
    -> std::pair<std::shared_ptr<SampleManager>, AssocInput>;

}  // namespace gelex

#endif  // GELEX_ESTIMATOR_FREQ_REML_H_
