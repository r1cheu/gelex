#pragma once
#include <armadillo>

namespace chenx
{
using arma::dmat;
using arma::sp_dmat;
bool CheckIdentity(const dmat& inputs);
bool CheckIdentity(const sp_dmat& inputs);
std::string ToLowercase(std::string_view input);
}  // namespace chenx
