#pragma once
#include "chenx/optim/ai_updater.h"
#include "chenx/optim/em_updater.h"
#include "chenx/optim/fisher_updater.h"
#include "chenx/optim/matrix_updater.h"
#include "chenx/optim/nr_updater.h"
#include "loop.h"

namespace chenx {
namespace optim {
using namespace arma;
template <typename eT>
RemlLoop<eT>::RemlLoop(const Col<eT>& y, const Mat<eT>& X, const Cube<eT>& Z, const Cube<eT>& rand) {
    _zkztr = cal_zkztr(Z, rand);
    _base_std = stddev(y);
    _base_var = _base_std * _base_std;
    _y = (y - mean(y)) / _base_std;
    _X = X;
    _var(Col<eT>(_zkztr.n_slices, fill::ones) * 0.1);
}

template <typename eT>
void RemlLoop<eT>::_init_var_updater(std::string_view method, std::unique_ptr<VarianceUpdater<eT>>& var_updater) {
    if (method == "ai") {
        var_updater = std::make_unique<AIUpdater<eT>>(_var, _y);
    } else if (method == "nr") {
        var_updater = std::make_unique<NRUpdater<eT>>(_var, _y);
    } else if (method == "fisher") {
        var_updater = std::make_unique<FisherUpdater<eT>>(_var, _y);
    } else {
        throw std::invalid_argument("Invalid method");
    }
}

template <typename eT>
eT RemlLoop<eT>::_cal_loglik(const Mat<eT>& v, const Mat<eT>& txvx, const Mat<eT>& proj_y) {
    return -0.5 * (log_det_sympd(v) + log_det_sympd(txvx) + as_scalar(_y.t() * proj_y));
}

template <typename eT>
bool RemlLoop<eT>::_has_converged(eT log_diff, eT tolerance) {
    if (std::abs(log_diff) < tolerance) {
        std::cout << "Converged!!!" << std::endl;
        return true;
    }
    return false;
}

template <typename eT>
void RemlLoop<eT>::run(std::string_view method, bool em_init, size_t max_iteration, double tolerance) {
    MatrixUpdater<eT> mat_updater(_X, _y, _zkztr);

    mat_updater.update(_var);

    eT logL = _cal_loglik(mat_updater.getV(), mat_updater.getTxvx(), mat_updater.getProj_y());
    eT log_diff = 0;

    size_t init_iteration = 1;
    Logger<eT> logger;
    if (em_init) {
        logger.start();
        // EM initialization of variance
        EMUpdater<eT> em_updater(_var, _y);
        _var = em_updater.update(mat_updater.getProj_y(), mat_updater.getPdv());

        mat_updater.update(_var);
        eT new_logL = _cal_loglik(mat_updater.getV(), mat_updater.getTxvx(), mat_updater.getProj_y());
        log_diff = new_logL - logL;

        logger.log(init_iteration, "EM", new_logL, _var);
        init_iteration++;
    }

    std::unique_ptr<VarianceUpdater<eT>> var_updater;
    _init_var_updater(method, var_updater);

    for (size_t i{init_iteration}; i < max_iteration; i++) {
        logger.start();

        _var = var_updater->update(mat_updater.getProj_y(), mat_updater.getPdv());
        mat_updater.update(_var);

        eT new_logL = _cal_loglik(mat_updater.getV(), mat_updater.getTxvx(), mat_updater.getProj_y());
        log_diff = new_logL - logL;
        logger.log(i, method, new_logL, _var);

        if (_has_converged(log_diff, tolerance)) {
            converged = true;
            break;
        }
    }

    if (!converged) {
        std::cout << "Not converged!!! Try to increase num of iteration" << std::endl;
    }
}

template <typename eT>
bool check_identity(const Mat<eT>& inputs) {
    if (inputs.n_rows != inputs.n_cols) {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i) {
        for (size_t j = 0; j < inputs.n_cols; ++j) {
            if (i != j && inputs(i, j) != 0) {
                return false;
            } else {
                if (i == j && inputs(i, j) != 1) {
                    return false;
                }
            }
        }
    }
    return true;
}

template <typename eT>
Mat<eT> cal_zkz(const Mat<eT>& z, const Col<eT>& k) {
    bool z_identity = check_identity(z);
    bool k_identity = check_identity(k);

    if (k_identity) {
        if (z_identity) {
            return k;
        } else {
            return z * z.t();
        }
    } else {
        if (z_identity) {
            return k.t() * k;
        } else {
            return z * k * z.t();
        }
    }
}

template <typename eT>
Cube<eT> cal_zkztr(const Cube<eT>& z, const Cube<eT>& k) {
    Cube<eT> result(z.n_rows, z.n_rows, k.n_slices + 1, fill::zeros);
    for (size_t i = 0; i < k.n_slices; ++i) {
        result.slice(i) = cal_zkz(z.slice(i), k.slice(i));
    }

    result.slice(k.n_slices + 1) = eye<Mat<eT>>(z.n_rows);
    return result;
}

}  // namespace optim
}  // namespace chenx
