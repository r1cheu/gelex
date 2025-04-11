#include <cstdint>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <armadillo>
#include <utility>

#include "gelex/python/dense_caster.h"

namespace bind
{
namespace nb = nanobind;
using namespace nb::literals;

class A
{
   public:
    explicit A(arma::dmat arr) : arr_(std::move(arr)) {}

    const arma::dmat& arr() const { return arr_; }
    arma::dmat& arr() { return arr_; }
    arr2d new_arr()
    {
        arma::dmat tmp(10, 10, arma::fill::randu);
        arr_ptr_ = reinterpret_cast<intptr_t>(tmp.memptr());
        return to_py(std::move(tmp));
    }
    intptr_t new_arr_ptr() const { return arr_ptr_; }

   private:
    arma::dmat arr_;
    intptr_t arr_ptr_;
};

intptr_t get_vec_ptr_value(arr1d arr)
{
    arma::dvec vec = to_arma(arr);
    return reinterpret_cast<intptr_t>(vec.memptr());
}

intptr_t get_mat_ptr_value(arr2d arr)
{
    arma::dmat mat = to_arma(arr);
    return reinterpret_cast<intptr_t>(mat.memptr());
}

std::optional<intptr_t> get_optional_mat_ptr_value(std::optional<arr2d> arr)
{
    if (!arr)
        return std::nullopt;
    arma::dmat mat = to_arma(*arr);
    return reinterpret_cast<intptr_t>(mat.memptr());
}

A return_A()
{
    A a(arma::dmat(3, 3, arma::fill::randu));
    a.arr()(0, 0) = 99.0;
    return a;
}

arr2d return_mat_rvalue()
{
    arma::dmat new_mat(3, 3, arma::fill::randu);
    new_mat(0, 0) = 99.0;
    return to_py(std::move(new_mat));
}

arr2d create_test_matrix(double value)
{
    arma::dmat mat(3, 3, arma::fill::ones);
    mat *= value;
    return to_py(std::move(mat));
}

arr2d modify_and_return_matrix(arr2d input)
{
    arma::dmat mat = to_arma(input);
    mat *= 2.0;
    return input;
}

void register_dense_caster_test(nb::module_& m)
{
    m.def("get_vec_ptr_value", &get_vec_ptr_value);
    m.def("get_mat_ptr_value", &get_mat_ptr_value);
    m.def(
        "get_optional_mat_ptr_value",
        &get_optional_mat_ptr_value,
        "arr"_a = nb::none());
    m.def("return_A", &return_A);
    m.def("return_mat_rvalue", &return_mat_rvalue);
    m.def("create_test_matrix", &create_test_matrix);
    m.def("modify_and_return_matrix", &modify_and_return_matrix);
}

NB_MODULE(_test, m)
{
    register_dense_caster_test(m);

    nb::class_<A>(m, "A")
        .def(
            "__init__",
            [](A* self, arr2d arr) { new (self) A(to_arma(std::move(arr))); })
        .def(
            "get_arr_ptr",
            [](const A& self)
            { return reinterpret_cast<intptr_t>(self.arr().memptr()); })
        .def("new_arr_ptr", &A::new_arr_ptr)
        .def("new_arr", &A::new_arr, nb::rv_policy::reference)
        .def(
            "const_arr",
            [](const A& self) { return to_py(self.arr()); },
            nb::rv_policy::reference_internal)
        .def(
            "mutable_arr",
            [](A& self) { return to_py(self.arr()); },
            nb::rv_policy::reference_internal);
}

}  // namespace bind
//
