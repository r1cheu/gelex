// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#ifndef GELEX_EXCEPTION_H_
#define GELEX_EXCEPTION_H_

#include <stdexcept>

namespace gelex
{

class GelexException : public std::runtime_error
{
   public:
    using std::runtime_error::runtime_error;
};

}  // namespace gelex

#endif  // GELEX_EXCEPTION_H_
