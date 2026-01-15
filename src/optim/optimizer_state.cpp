#include "gelex/optim/optimizer_state.h"

#include "gelex/model/freq/model_new.h"

namespace gelex
{

OptimizerState::OptimizerState(const FreqModel& model)
    : num_individuals_(model.num_individuals()),
      phenotype_var_(
          model.phenotype().squaredNorm()
          / static_cast<double>(num_individuals_ - 1))
{
    v.resize(num_individuals_, num_individuals_);
    proj.resize(num_individuals_, num_individuals_);
    proj_y.resize(num_individuals_);
    tx_vinv_x.resize(model.fixed().X.cols(), model.fixed().X.cols());
}

}  // namespace gelex
