#pragma once
#include "variance_updater.h"
namespace chenx {
namespace optim {
template <typename eT>
class FisherUpdater : public VarianceUpdater<eT> {
   public:
    FisherUpdater(Col<eT> init_var, const Col<eT>& y) : VarianceUpdater<eT>(init_var, y) {}

   private:
    virtual eT _cal_info_element(const Mat<eT>& proj_y, const Mat<eT>& pdvi, const Mat<eT>& pdvj) override {
        return -0.5 * trace(pdvi * pdvj);
    }
};
}  // namespace optim
}  // namespace chenx
