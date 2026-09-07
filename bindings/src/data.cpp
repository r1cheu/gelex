// Copyright 2026 RuLei Chen
// SPDX-License-Identifier: Apache-2.0

#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>

#include "gelex/data/encode/matrix.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

void register_data(nb::module_& m)
{
    m.def(
        "encode_inplace",
        &gelex::encode_inplace,
        nb::arg("genotypes").noconvert(),
        nb::arg("mode"),
        nb::arg("method"));
}

}  // namespace gelex_py
