#pragma once
#include "variance_updater.h"
namespace chenx {
namespace optim {
template <typename eT>
class NRUpdater : public VarianceUpdater<eT> {
   public:
    NRUpdater(Col<eT> init_var, const Col<eT>& y) : VarianceUpdater<eT>(init_var, y) {}

   private:
    virtual eT _cal_info_element(const Mat<eT>& proj_y, const Mat<eT>& pdvi, const Mat<eT>& pdvj) override {
        return 0.5 * trace(pdvi * pdvj) - as_scalar(this->_y.t() * pdvi * pdvj * proj_y);
    }
};
}  // namespace optim
}  // namespace chenx
