#pragma once
#include "armadillo"
#include "variance_updater.h"

namespace chenx
{

using namespace arma;
template <typename eT>
VarianceUpdater<eT>::VarianceUpdater(Col<eT>& init_var, const Col<eT>& y) : _var{init_var}, _y{y}
{
    auto n = _var.n_elem;
    y_var = var(y);
    _prev_var.zeros(n);
    _score.zeros(n);
    _info_matrix.zeros(n, n);
    _info_matrix_inv.zeros(n, n);
}

template <typename eT>
void VarianceUpdater<eT>::_constrain_var(Col<eT>& var)
{
    double delta = 0, constr_scale = 1e-6;
    Col<eT> constrained(var.n_elem, fill::zeros);

    for (size_t i = 0; i < var.n_elem; i++)
    {
        if (var(i) < 0)
        {
            delta += y_var * constr_scale - var(i);
            var(i) = y_var * constr_scale;
            constrained(i) = 1;
        }
    }
    delta /= var.n_elem - sum(constrained);

    for (size_t i = 0; i < var.n_elem; i++)
    {
        if (constrained(i) == 0 && var(i) > delta)
        {
            var(i) -= delta;
        }
    }

    if (sum(constrained) > var.n_elem / 2)
    {
        std::cout << "half of the variance components are constrained! The estimate is not reliable. " << std::endl;
    }
}

template <typename eT>
void VarianceUpdater<eT>::_cal_score(const Mat<eT>& proj_y, const Cube<eT>& pdv)
{
    _score(0) = -0.5 * (trace(pdv.slice(0)) - as_scalar(_y.t() * pdv.slice(0) * proj_y));
    for (size_t i{1}; i < _var.n_elem; ++i)
    {
        _score(i) = -0.5 * (trace(pdv.slice(i)) - as_scalar(_y.t() * pdv.slice(i) * proj_y));
    }
}

template <typename eT>
void VarianceUpdater<eT>::_cal_info_matrix(const Mat<eT>& proj_y, const Cube<eT>& pdv)
{
    auto n = pdv.n_slices;
    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            _info_matrix(i, j) = _cal_info_element(proj_y, pdv.slice(i), pdv.slice(j));
            if (i != j)
            {
                _info_matrix(j, i) = _info_matrix(i, j);
            }
        }
    }
}

template <typename eT>
Col<eT> VarianceUpdater<eT>::update(const Mat<eT>& proj_y, const Cube<eT>& pdv, double lambda)
{
    _prev_var = _var;
    _cal_score(proj_y, pdv);
    _cal_info_matrix(proj_y, pdv);
    if (!pinv(_info_matrix_inv, _info_matrix))
    {
        throw std::runtime_error("Matrix inversion failed");
    }
    Col<eT> delta = -lambda * _info_matrix_inv * _score;
    _var += delta;
    _constrain_var(_var);

    return Col<eT>(_var);
}

template <typename eT>
eT VarianceUpdater<eT>::get_vardiff()
{
    return norm(_var - _prev_var) / norm(_var);
}

}  // namespace chenx
