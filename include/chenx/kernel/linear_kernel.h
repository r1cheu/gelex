#pragma once
#include "armadillo"
#include "kernel.h"

namespace chenx {

using namespace arma;

template <typename eT>
class LinearKernel : public Kernel<eT> {
   protected:
    eT compute(const Col<eT>& i, const Col<eT>& j) const override { return dot(i, j); }
};

}  // namespace chenx
