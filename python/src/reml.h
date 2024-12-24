#pragma once
// clang-format off
#include "carma"
// clang-format on
#include <chenx/optim/variance_updater.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace binds
{
using namespace arma;
namespace py = pybind11;
using namespace chenx;

template <typename eT>
class REMLLoop
{
  public:
    REMLLoop(
        const py::array_t<eT>& y_arr,
        const py::array_t<eT>& X_arr,
        const py::array_t<uword>& z_index_arr,
        const py::array_t<eT>& rands_arr,
        std::vector<std::string> rand_names);
    void run(
        const py::array_t<eT>& varcomp,
        std::string_view method,
        bool em_init,
        size_t max_iteration,
        double tolerance,
        bool verbose);
    py::array_t<eT> get_varcomp() const;
    py::array_t<eT> get_beta() const;
    py::array_t<eT> get_blup() const;
    py::array_t<eT> get_gebv(const py::array_t<eT>& full_X) const;

  private:
    Cube<eT> _zkztr;
    const Col<eT> _y;
    const Mat<eT> _X;
    const Cube<eT> _rands;
    std::vector<SpMat<eT>> _Z;
    std::vector<std::string> _rand_names;
    bool _converged = false;
    void _init_varcomp(Col<eT>& varcomp);
    void _init_var_updater(
        Col<eT>& varcomp,
        std::string_view method,
        std::unique_ptr<VarianceUpdater<eT>>& var_updater);
    eT _cal_loglik(
        const double& logdet_v,
        const Mat<eT>& txvx,
        const Mat<eT>& proj_y);
    bool _has_converged(eT var_diff, eT log_diff, eT tolerance);
    void set_blup(const Mat<eT>& Vi);
    void set_beta(const Mat<eT>& txvx, const Mat<eT>& Vi);
    Col<eT> _varcomp;
    Col<eT> _beta;
    Mat<eT> _blup;
};
} // namespace binds

#include "reml_impl.h"
