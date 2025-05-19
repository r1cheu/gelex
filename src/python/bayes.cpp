#include <optional>

#include <fmt/format.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
#include <armadillo>

#include "gelex/model/bayes_model.h"
#include "gelex/model/bayes_prior.h"
#include "gelex/python/dense.h"
#include "gelex/python/sparse.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
namespace gx = gelex;

using arma::dmat;
using arma::dvec;
template <typename BayesModel>
void register_bayes_model(nb::module_& module, const char* name)
{
    nb::class_<BayesModel>(module, name)
        .def(
            nb::init<dvec, dmat, std::optional<dmat>>(),
            "phenotype"_a.noconvert(),
            "genotype_mat"_a.noconvert(),
            "design_mat_beta"_a.noconvert() = nb::none(),
            nb::keep_alive<1, 2>(),
            nb::keep_alive<1, 3>(),
            nb::keep_alive<1, 4>())
        .def(
            "priors",
            [](BayesModel& self) { return self.priors(); },
            nb::rv_policy::reference_internal)
        .def(
            "__repr__",
            [](const BayesModel& self)
            {
                return fmt::format(
                    "<{} object at {:p}: phenotype({:d}, 1), "
                    "genotype_mat({:d}, {:d}), "
                    "design_mat_beta({:d}, {:d})",
                    BayesModel::name,
                    static_cast<const void*>(&self),
                    self.phenotype().n_elem,
                    self.genotype_mat().n_rows,
                    self.genotype_mat().n_cols,
                    self.design_mat_beta() ? self.design_mat_beta()->n_rows : 0,
                    self.design_mat_beta() ? self.design_mat_beta()->n_cols
                                           : 0);
            })
        .def(
            "__str__",
            [](const BayesModel& self)
            {
                std::string info = fmt::format(
                    "┌─ {} Model ─────────────────────────────────\n"
                    "│ Individuals:        {:6d}\n",
                    BayesModel::name,
                    self.phenotype().n_elem);

                // Show fixed effects only if they exist
                if (self.design_mat_beta())
                {
                    info += fmt::format(
                        "│ Fixed Effects:      {:6d}\n",
                        self.design_mat_beta()->n_cols);
                }
                // Show environmental effects only if they exist
                for (size_t i = 0; i < self.design_mat_r().size(); ++i)
                {
                    info += fmt::format(
                        "│ Random Eff. {}:       {:6d}\n",
                        self.random_names()[i],
                        self.design_mat_r()[i].n_cols);
                }
                info += fmt::format(
                    "│ SNPs:               {:6d}\n",
                    self.genotype_mat().n_cols);

                info += "├─ Priors ──────────────────────────────────────\n";
                info += fmt::format(
                    "│ Variance Priors:\n"
                    "│   σₐ (additive):      nu = {:6.1f}, s² = {:6.4f}\n",
                    self.priors().sigma_a.nu,
                    self.priors().sigma_a.s2);

                if (self.design_mat_r().size() > 0)
                {
                    info += fmt::format(
                        "│   σᵣ (random): nu = {:6.1f}, s² = {:6.4f}\n",
                        self.priors().sigma_r.nu,
                        self.priors().sigma_r.s2);
                }

                info += fmt::format(
                    "│   σₑ (residual):      nu = {:6.1f}, s² = {:6.4f}\n",
                    self.priors().sigma_e.nu,
                    self.priors().sigma_e.s2);

                if (BayesModel::has_pi)
                {
                    info += "│ Mixture Priors";
                    if (!BayesModel::fixed_pi)
                    {
                        info += " (Fixed)";
                    }
                    info += ":\n";
                    info += fmt::format(
                        "│   π = [{}]", fmt::join(self.priors().pi, ", "));
                    info += "\n";
                }
                info += "└───────────────────────────────────────────────";
                return info;
            });
}

void bayesa(nb::module_& module)
{
    register_bayes_model<gx::BayesA>(module, "BayesA");
}

void bayesrr(nb::module_& module)
{
    register_bayes_model<gx::BayesRR>(module, "BayesRR");
}

void bayesb(nb::module_& module)
{
    register_bayes_model<gx::BayesB>(module, "BayesB");
}

void bayesbpi(nb::module_& module)
{
    register_bayes_model<gx::BayesBpi>(module, "BayesBpi");
}

void bayesc(nb::module_& module)
{
    register_bayes_model<gx::BayesC>(module, "BayesC");
}

void bayescpi(nb::module_& module)
{
    register_bayes_model<gx::BayesCpi>(module, "BayesCpi");
}

void sigma_prior(nb::module_& module)
{
    nb::class_<gx::sigma_prior>(module, "sigma_prior")
        .def(nb::init<double, double>(), "nu"_a, "s2"_a)
        .def_rw("nu", &gx::sigma_prior::nu, nb::rv_policy::reference_internal)
        .def_rw("s2", &gx::sigma_prior::s2, nb::rv_policy::reference_internal);
}

void priors(nb::module_& module)
{
    nb::class_<gelex::Priors>(module, "Priors")
        .def(nb::init<>())
        .def_rw("pi", &gx::Priors::pi, nb::rv_policy::reference_internal)
        .def_rw(
            "sigma_a", &gx::Priors::sigma_a, nb::rv_policy::reference_internal)
        .def_rw(
            "sigma_r", &gx::Priors::sigma_r, nb::rv_policy::reference_internal)
        .def_rw(
            "sigma_e", &gx::Priors::sigma_e, nb::rv_policy::reference_internal);
}

}  // namespace bind
