#pragma once
#include <armadillo>
#include <memory>
#include <string_view>
#include "chenx/optim/variance_updater.h"
#include "log.h"

namespace chenx {
namespace optim {
using namespace arma;
template <typename eT>
class RemlLoop {
   public:
    RemlLoop(const Col<eT>& y, const Mat<eT>& X, const Cube<eT>& Z, const Cube<eT>& rand);
    void run(std::string_view method, bool em_init, size_t max_iteration, double tolerance);

   private:
    const Cube<eT> _zkztr;
    const Col<eT> _y;
    const Mat<eT> _X;
    Col<eT> _var;
    double _base_std;
    double _base_var;
    bool converged = false;
    void _init_var_updater(std::string_view method, std::unique_ptr<VarianceUpdater<eT>>& var_updater);
    eT _cal_loglik(const Mat<eT>& v, const Mat<eT>& txvt, const Mat<eT>& proj_y);
    bool _has_converged(eT log_diff, eT tolerance);
};

template <typename eT>
bool check_identity(const Mat<eT>& inputs);

template <typename eT>
Mat<eT> cal_zkz(const Mat<eT>& z, const Col<eT>& k);

template <typename eT>
Cube<eT> cal_zkztr(const Cube<eT>& z, const Cube<eT>& k);
}  // namespace optim
}  // namespace chenx

#include "loop_impl.h"
