#pragma once
#include <armadillo>

namespace chenx {
namespace optim {
using namespace arma;
template <typename eT>
class VarianceUpdater {
   public:
    VarianceUpdater(Col<eT> init_var, const Col<eT>& y);
    virtual ~VarianceUpdater() = default;
    Col<eT> update(const Mat<eT>& proj_y, const Cube<eT>& pdv, double lambda = 1.0);

   protected:
    const Col<eT>& _y;

   private:
    Col<eT> _var;
    Col<eT> _score;
    Mat<eT> _info_matrix;
    Mat<eT> _info_matrix_inv;
    void _cal_score(const Mat<eT>& proj_y, const Cube<eT>& pdv);
    void _cal_info_matrix(const Mat<eT>& proj_y, const Cube<eT>& pdv);
    virtual eT _cal_info_element(const Mat<eT>& proj_y, const Mat<eT>& pdvi, const Mat<eT>& pdvj) = 0;
};
}  // namespace optim
}  // namespace chenx

#include "variance_updater_impl.h"
