#include "chenx/utils.h"
#include <cmath>

namespace chenx
{
bool CheckIdentity(const dmat& inputs)
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

bool CheckIdentity(const sp_dmat& inputs)
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

std::string ToLowercase(std::string_view input)
{
    std::string result(input.begin(), input.end());
    std::transform(result.begin(), result.end(), result.begin(), ::tolower);
    return result;
}

}  // namespace chenx
