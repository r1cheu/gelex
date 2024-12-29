#pragma once
#include <armadillo>
#include <cmath>

namespace chenx
{
using namespace arma;
bool check_identity(const dmat& inputs);
bool check_identity(const sp_dmat& inputs);
}  // namespace chenx
