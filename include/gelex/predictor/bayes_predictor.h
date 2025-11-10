#pragma once
#include <expected>
#include <string_view>

#include "gelex/error.h"

class BayesPredictor
{
   public:
    static auto create(std::string_view bed, std::string_view result_prefix);

   private:
};
