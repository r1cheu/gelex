#include "gelex/python/dense.h"
#include "gelex/python/sparse.h"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <armadillo>
#include <optional>

namespace nb = nanobind;
using namespace nb::literals;
using namespace arma;

namespace bind
{
namespace nb = nanobind;
using namespace nb::literals;

struct ClassWithArmaMember
{
    dmat member = dmat(2, 2, fill::ones);
    const dmat& get_member_ref() { return member; }
    const dmat get_member_copy() { return member; }
};

class ClassInitFromPython
{
   public:
    ClassInitFromPython(dmat mat) : mat_{std::move(mat)} {}
    dmat& get_mat() { return mat_; }

   private:
    dmat mat_;
};

class ClassInitOptional
{
   public:
    ClassInitOptional(std::optional<dmat> mat) : mat_{std::move(mat)} {}
    dmat& get_mat()
    {
        if (mat_)
        {
            return *mat_;
        }
        return mat_default_;
    }

   private:
    std::optional<dmat> mat_;
    dmat mat_default_{{0, 2, 3}, {1, 2, 3}};
};

void register_dense_test(nb::module_& m)
{
    m.def(
        "add_ivec",
        [](const ivec& a, const ivec& b) -> ivec { return a + b; },
        "a"_a,
        "b"_a.noconvert());
    m.def(
        "add_imat",
        [](const imat& a, const imat& b) -> imat { return a + b; },
        "a"_a,
        "b"_a.noconvert());
    m.def("update_ivec", [](ivec& a) { a(0) = 99.0; }, "a"_a.noconvert());
    m.def("update_imat", [](imat& a) { a(0, 0) = 99; }, "a"_a.noconvert());

    nb::class_<ClassWithArmaMember>(m, "ClassWithArmaMember")
        .def(nb::init<>())
        .def_prop_ro("member_ro_ref", &ClassWithArmaMember::get_member_ref)
        .def_prop_ro("member_ro_copy", &ClassWithArmaMember::get_member_copy)
        .def_rw("member", &ClassWithArmaMember::member);

    nb::class_<ClassInitFromPython>(m, "ClassInitFromPython")
        .def(nb::init<dmat>(), nb::keep_alive<1, 2>(), "mat"_a.noconvert())
        .def(
            "get_mat",
            &ClassInitFromPython::get_mat,
            nb::rv_policy::reference_internal);

    nb::class_<ClassInitOptional>(m, "ClassInitOptional")
        .def(
            nb::init<std::optional<dmat>>(),
            nb::keep_alive<1, 2>(),
            "mat"_a.noconvert() = nb::none())
        .def(
            "get_mat",
            &ClassInitOptional::get_mat,
            nb::rv_policy::reference_internal);
    m.def(
        "sparse",
        []() { return sp_dmat(dmat{{1, 3, 4}, {0, 0, 0}, {1, 0, 0}}); },
        nb::rv_policy::move);
    m.def(
        "sparse_add",
        [](sp_dmat a) -> dmat { return a + 1; },
        nb::rv_policy::move);
}
}  // namespace bind
NB_MODULE(_test, m)
{
    bind::register_dense_test(m);
}
