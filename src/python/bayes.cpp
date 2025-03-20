#include <optional>

#include <fmt/format.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <armadillo>

#include "gelex/model/bayes_model.h"
#include "gelex/python/array_caster.h"
namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;

template <typename BayesModel>
void register_bayes_model(nb::module_& module, const char* name)
{
    nb::class_<BayesModel>(module, name)
        .def(
            "__init__",
            [](BayesModel* self,
               arr1d&& phenotype,
               arr2d&& genotype_mat,
               std::optional<arr2d> design_mat_beta,
               std::optional<arr2d> design_mat_r)
            {
                new (self) BayesModel{
                    ToArma(std::move(phenotype)),
                    ToArma(std::move(genotype_mat)),
                    design_mat_beta
                        ? std::make_optional(ToArma(*design_mat_beta))
                        : std::nullopt,
                    design_mat_r ? std::make_optional(ToArma(*design_mat_r))
                                 : std::nullopt};
            },
            "phenotype"_a.noconvert(),
            "genotype_mat"_a.noconvert(),
            "design_mat_beta"_a.noconvert() = nb::none(),
            "design_mat_r"_a.noconvert() = nb::none(),
            nb::keep_alive<1, 2>(),
            nb::keep_alive<1, 3>(),
            nb::keep_alive<1, 4>(),
            nb::keep_alive<1, 5>())
        .def(
            "__repr__",
            [](const BayesModel& self)
            {
                return fmt::format(
                    "<{} object at {:p}: phenotype({:d}, 1), "
                    "genotype_mat({:d}, {:d}), "
                    "design_mat_beta({:d}, {:d}), "
                    "design_mat_r({:d}, {:d})>",
                    BayesModel::name,
                    static_cast<const void*>(&self),
                    self.phenotype().n_elem,
                    self.genotype_mat().n_rows,
                    self.genotype_mat().n_cols,
                    self.design_mat_beta() ? self.design_mat_beta()->n_rows : 0,
                    self.design_mat_beta() ? self.design_mat_beta()->n_cols : 0,
                    self.design_mat_r() ? self.design_mat_r()->n_rows : 0,
                    self.design_mat_r() ? self.design_mat_r()->n_cols : 0);
            })
        .def(
            "__str__",
            [](const BayesModel& self)
            {
                return fmt::format(
                    "{}\n{:d} Individuals, {:d} Fixed Effects, "
                    "{:d} Environmental Effects, {:d} SNPs",
                    BayesModel::name,
                    self.phenotype().n_elem,
                    self.design_mat_beta() ? self.design_mat_beta()->n_cols : 0,
                    self.design_mat_r() ? self.design_mat_r()->n_cols : 0,
                    self.genotype_mat().n_cols);
            });
}

void bayesa(nb::module_& module)
{
    register_bayes_model<gelex::BayesA>(module, "BayesA");
}

void bayesrr(nb::module_& module)
{
    register_bayes_model<gelex::BayesRR>(module, "BayesRR");
}

void bayesb(nb::module_& module)
{
    register_bayes_model<gelex::BayesB>(module, "BayesB");
}

void bayesbpi(nb::module_& module)
{
    register_bayes_model<gelex::BayesBpi>(module, "BayesBpi");
}

void bayesc(nb::module_& module)
{
    register_bayes_model<gelex::BayesC>(module, "BayesC");
}

void bayescpi(nb::module_& module)
{
    register_bayes_model<gelex::BayesCpi>(module, "BayesCpi");
}

}  // namespace bind
