#pragma once
#include <armadillo>

namespace chenx
{
using namespace arma;
template <typename eT>
class MatrixUpdater
{
  public:
    MatrixUpdater(const Mat<eT>& X, const Col<eT>& y, const Cube<eT>& zkztr);
    void update(const Col<eT>& var);
    const Mat<eT>& getVi() const
    {
        return _v;
    }
    const eT getLogdet_v() const
    {
        return _log_det_v;
    }
    const Col<eT>& getProj_y() const
    {
        return _proj_y;
    }
    const Cube<eT>& getPdv() const
    {
        return _pdv;
    }
    const Mat<eT>& getTxvx() const
    {
        return _txvx;
    }

    static eT inv_log_det_sympd(Mat<eT>& V);

  private:
    const Mat<eT>& _X;
    const Col<eT>& _y;
    const Cube<eT>& _zkztr;
    double _log_det_v;
    double _inv_tol;
    Col<eT> _proj_y;
    Mat<eT> _v, _proj, _txvx;
    Cube<eT> _pdv;
    void _cal_v(const Col<eT>& var);
    void _cal_proj_matrix();
    void _cal_pdv();
};

} // namespace chenx

#include "matrix_updater_impl.h"
