#include <nanobind/nanobind.h>

#include "gelex_py/register.h"

namespace nb = nanobind;

NB_MODULE(_gelex, m)
{
    m.doc() = "Python bindings for the gelex C++ library";

    gelex_py::register_types(m);
}
