#pragma once

// clang-format off
#include "reml.h"
// clang-format on
#include <chenx/log/log.h>
#include <chenx/optim/em_updater.h>
#include <chenx/optim/gradient_calculater.h>
#include <chenx/optim/variance_updater.h>
#include <chenx/optim/zkztr.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace binds
{
using namespace arma;
namespace py = pybind11;
using namespace chenx;
template <typename eT>
REMLLoop<eT>::REMLLoop(
    const py::array_t<eT>& y_arr,
    const py::array_t<eT>& X_arr,
    const py::array_t<uword>& z_index_arr,
    const py::array_t<eT>& rands_arr,
    std::vector<std::string> rand_names)
    : _y{carma::arr_to_col_view(y_arr)},
      _X{carma::arr_to_mat_view(X_arr)},
      _rands{carma::arr_to_cube_view(rands_arr)},
      _rand_names{rand_names}
{
    Col<uword> z_index = carma::arr_to_col(z_index_arr);
    _Z = create_z<eT>(_rands.n_slices, z_index, _rands.n_rows);
    _zkztr = cal_zkztr(_Z, _rands);
    _varcomp = Col<eT>(_zkztr.n_slices);
    _beta = Col<eT>(_X.n_cols);
    _blup = Mat<eT>(
        _Z[0].n_cols, static_cast<uword>(_Z.size()), fill::value(datum::nan));
}

template <typename eT>
void REMLLoop<eT>::_init_var_updater(
    Col<eT>& varcomp,
    std::string_view method,
    std::unique_ptr<VarianceUpdater<eT>>& var_updater)
{
    if (method == "AI")
    {
        var_updater = std::make_unique<AIUpdater<eT>>(varcomp, _y);
    }
    else if (method == "NR")
    {
        var_updater = std::make_unique<NRUpdater<eT>>(varcomp, _y);
    }
    else if (method == "Fisher")
    {
        var_updater = std::make_unique<FisherUpdater<eT>>(varcomp, _y);
    }
    else
    {
        throw std::invalid_argument(
            "Invalid optim method, please choose from `AI`, `NR`, `Fisher`");
    }
}

template <typename eT>
void REMLLoop<eT>::_init_varcomp(Col<eT>& varcomp)
{
    varcomp /= sum(varcomp);
    varcomp *= as_scalar(cov(_y));
}

template <typename eT>
eT REMLLoop<eT>::_cal_loglik(
    const double& logdet_v,
    const Mat<eT>& txvx,
    const Mat<eT>& proj_y)
{
    return -0.5 * (logdet_v + log_det_sympd(txvx) + as_scalar(_y.t() * proj_y));
}

template <typename eT>
bool REMLLoop<eT>::_has_converged(eT var_diff, eT log_diff, eT tolerance)
{
    if (var_diff < tolerance)
    {
        eT abs_log_diff = fabs(log_diff);
        if (abs_log_diff < 1e-4 || (abs_log_diff < 1e-2 && log_diff < 0))
        {
            return true;
        }
    }
    return false;
}

template <typename eT>
void REMLLoop<eT>::run(
    const py::array_t<eT>& varcomp_arr,
    std::string_view method,
    bool em_init,
    size_t max_iteration,
    double tolerance,
    bool verbose)
{
    Col<eT> varcomp = carma::arr_to_col(varcomp_arr);
    MatrixUpdater<eT> mat_updater(_X, _y, _zkztr);
    _init_varcomp(varcomp);
    mat_updater.update(varcomp);
    eT logL = _cal_loglik(
        mat_updater.getLogdet_v(),
        mat_updater.getTxvx(),
        mat_updater.getProj_y());
    eT log_diff = 1e100;
    size_t iteration = 1;
    Logger<eT> logger(_rand_names, verbose);

    if (em_init)
    {
        logger.start();
        EMUpdater<eT> em_updater(varcomp, _y);
        varcomp
            = em_updater.update(mat_updater.getProj_y(), mat_updater.getPdv());
        mat_updater.update(varcomp);
        logL = _cal_loglik(
            mat_updater.getLogdet_v(),
            mat_updater.getTxvx(),
            mat_updater.getProj_y());
        logger.log(iteration++, "EM", logL, varcomp);
    }

    std::unique_ptr<VarianceUpdater<eT>> var_updater;
    _init_var_updater(varcomp, method, var_updater);
    for (; iteration < max_iteration; ++iteration)
    {
        logger.start();

        varcomp = var_updater->update(
            mat_updater.getProj_y(), mat_updater.getPdv());
        mat_updater.update(varcomp);

        eT new_logL = _cal_loglik(
            mat_updater.getLogdet_v(),
            mat_updater.getTxvx(),
            mat_updater.getProj_y());
        log_diff = new_logL - logL;
        logL = new_logL;
        logger.log(iteration, method, logL, varcomp);

        if (_has_converged(var_updater->get_vardiff(), log_diff, tolerance))
        {
            _converged = true;
            break;
        }
    }
    _varcomp = varcomp;
    logger.end();
    set_beta(mat_updater.getTxvx(), mat_updater.getVi());
    set_blup(mat_updater.getVi());
    std::cout
        << (_converged ? "Converged!!!"
                       : "Not converged!!! Try to increase num of iteration")
        << std::endl;
}

template <typename eT>
py::array_t<eT> REMLLoop<eT>::get_varcomp() const
{
    return carma::col_to_arr(_varcomp);
}

template <typename eT>
py::array_t<eT> REMLLoop<eT>::get_beta() const
{
    return carma::col_to_arr(_beta);
}

template <typename eT>
py::array_t<eT> REMLLoop<eT>::get_blup() const
{
    return carma::mat_to_arr(_blup);
}

template <typename eT>
void REMLLoop<eT>::set_blup(const Mat<eT>& Vi)
{
    Col<eT> fitted = _X * _beta;
    Col<eT> res = _y - fitted;
    Col<eT> vi_res = Vi * res;

    for (uword i = 0; i < _Z.size(); ++i)
    {
        _blup.col(i) = _rands.slice(i) * _varcomp(i) * _Z[i].t() * vi_res;
    }
}

template <typename eT>
void REMLLoop<eT>::set_beta(const Mat<eT>& txvx, const Mat<eT>& Vi)
{
    _beta = inv_sympd(txvx) * (_X.t() * Vi * _y);
}

template <typename eT>
py::array_t<eT> REMLLoop<eT>::get_gebv(const py::array_t<eT>& full_X_arr) const
{
    Mat<eT> full_X = carma::arr_to_mat_view(full_X_arr);
    return carma::col_to_arr(Col<eT>(full_X * _beta + sum(_blup, 1)));
}
}  // namespace binds
