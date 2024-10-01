#pragma once
#include <armadillo>

namespace chenx {
namespace optim {
using namespace arma;
template <typename eT>
class MatrixUpdater {
   public:
    MatrixUpdater(const Mat<eT>& X, const Col<eT>& y, const Cube<eT>& zkztr);
    void update(const Col<eT>& var);
    const Mat<eT>& getV() const { return _v; }
    const Mat<eT>& getVi() const { return _vi; }
    const Mat<eT>& getProj_y() const { return _proj_y; }
    const Mat<eT>& getPdv() const { return _pdv; }
    const Mat<eT>& getTxvx() const { return _txvx; }

   private:
    const Mat<eT>& _X;
    const Col<eT>& _y;
    const Cube<eT>& _zkztr;
    Col<eT> _proj_y;
    Mat<eT> _v, _vi, _proj, _txvx;
    Cube<eT> _pdv;
    void _cal_v(const Col<eT>& var);
    void _cal_proj_matrix();
    void _cal_pdv();
};
}  // namespace optim
}  // namespace chenx

#include "matrix_updater_impl.h"
