#pragma once
#include <armadillo>
#include "fmt/color.h"

namespace chenx
{
using arma::dmat;
using arma::dvec;
using arma::sp_dmat;
bool CheckIdentity(const dmat& inputs);
bool CheckIdentity(const sp_dmat& inputs);
std::string ToLowercase(std::string_view input);

template <typename T>
auto cyan(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::cyan));
}

template <typename T>
auto rebecca_purple(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::blue_violet));
}

template <typename T>
auto gray(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gray));
}

template <typename T>
auto gold(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::gold));
}

template <typename T>
auto red(const T& value)
{
    return fmt::styled(value, fmt::fg(fmt::color::red));
}

}  // namespace chenx
