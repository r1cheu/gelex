#pragma once
#include "chenx/utils.h"
#include <cmath>

namespace chenx
{
using namespace arma;
bool check_identity(const dmat& inputs)
{
    if (!inputs.is_square())
    {
        return false;
    }

    if (!inputs.is_diagmat())
    {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i)
    {
        if (inputs.at(i, i) != 1)
        {
            return false;
        }
    }
    return true;
}

bool check_identity(const sp_dmat& inputs)
{
    if (!inputs.is_square())
    {
        return false;
    }

    if (!inputs.is_diagmat())
    {
        return false;
    }

    for (size_t i = 0; i < inputs.n_rows; ++i)
    {
        if (inputs.at(i, i) != 1)
        {
            return false;
        }
    }
    return true;
}
}  // namespace chenx
