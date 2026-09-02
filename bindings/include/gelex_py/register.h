#ifndef GELEX_PY_REGISTER_H_
#define GELEX_PY_REGISTER_H_

#include <nanobind/nanobind.h>

namespace gelex_py
{

void register_types(nanobind::module_& m);
void register_data(nanobind::module_& m);

}  // namespace gelex_py

#endif  // GELEX_PY_REGISTER_H_
