#pragma once
#include "armadillo"
#include <utility>

namespace chenx
{
using namespace arma;
template <typename eT>
class EMUpdater
{
  public:
    EMUpdater(Col<eT>&& init_var, const Col<eT>& y);
    Col<eT> update(const Mat<eT>& proj_y, const Cube<eT>& pdv);

  private:
    const Col<eT>& _y;
    Col<eT> _var;
    double _N;
};

template <typename eT>
EMUpdater<eT>::EMUpdater(Col<eT>&& init_var, const Col<eT>& y)
    : _y{y}, _var{std::move(init_var)}
{
    _N = y.n_elem;
}

template <typename eT>
Col<eT> EMUpdater<eT>::update(const Mat<eT>& proj_y, const Cube<eT>& pdv)
{
    Col<eT> _var_2 = _var % _var;
    for (size_t i = 0; i < _var.n_elem; i++)
    {
        _var(i)
            = as_scalar(
                  _var_2.at(i) * (_y.t() * pdv.slice(i) * proj_y)
                  + trace(
                      _var.at(i) * eye(_N, _N) - _var_2.at(i) * pdv.slice(i)))
              / _N;
    }
    return Col<eT>(_var);
}
} // namespace chenx
