#pragma once
#include <armadillo>

namespace chenx
{
using namespace arma;
template <typename eT>
class VarianceUpdater
{
  public:
    VarianceUpdater(Col<eT>&& init_var, const Col<eT>& y);
    virtual ~VarianceUpdater() = default;
    Col<eT> update(
        const Mat<eT>& proj_y,
        const Cube<eT>& pdv,
        double lambda = 1.0);
    eT get_vardiff();
    Col<eT> get_var()
    {
        return _var;
    }
    const Col<eT>& get_y() const
    {
        return _y;
    }

  private:
    const Col<eT>& _y;
    Col<eT>& _var;
    double y_var;
    Col<eT> _prev_var;
    Col<eT> _score;
    Mat<eT> _info_matrix;
    Mat<eT> _info_matrix_inv;
    void _constrain_var(Col<eT>& var);
    void _cal_score(const Mat<eT>& proj_y, const Cube<eT>& pdv);
    void _cal_info_matrix(const Mat<eT>& proj_y, const Cube<eT>& pdv);
    virtual eT _cal_info_element(
        const Mat<eT>& proj_y,
        const Mat<eT>& pdvi,
        const Mat<eT>& pdvj)
        = 0;
};

// Fisher scoring algorithm
template <typename eT>
class FisherUpdater : public VarianceUpdater<eT>
{
  public:
    FisherUpdater(Col<eT>& init_var, const Col<eT>& y)
        : VarianceUpdater<eT>{init_var, y}
    {
    }

  private:
    virtual eT _cal_info_element(
        const Mat<eT>& proj_y,
        const Mat<eT>& pdvi,
        const Mat<eT>& pdvj) override
    {
        return -0.5 * trace(pdvi * pdvj);
    }
};

// Newton-Raphson algorithm
template <typename eT>
class NRUpdater : public VarianceUpdater<eT>
{
  public:
    NRUpdater(Col<eT>& init_var, const Col<eT>& y)
        : VarianceUpdater<eT>{init_var, y}
    {
    }

  private:
    virtual eT _cal_info_element(
        const Mat<eT>& proj_y,
        const Mat<eT>& pdvi,
        const Mat<eT>& pdvj) override
    {
        return 0.5 * trace(pdvi * pdvj)
               - as_scalar(
                   VarianceUpdater<eT>::get_y().t() * pdvi * pdvj * proj_y);
    }
};

// Average Information algorithm
template <typename eT>
class AIUpdater : public VarianceUpdater<eT>
{
  public:
    AIUpdater(Col<eT>& init_var, const Col<eT>& y)
        : VarianceUpdater<eT>{init_var, y}
    {
    }

  private:
    virtual eT _cal_info_element(
        const Mat<eT>& proj_y,
        const Mat<eT>& pdvi,
        const Mat<eT>& pdvj) override
    {
        return -0.5
               * as_scalar(
                   VarianceUpdater<eT>::get_y().t() * pdvi * pdvj * proj_y);
    }
};
} // namespace chenx

#include "variance_updater_impl.h"
