#include <string>
#include <string_view>
#include <vector>

#include <fmt/color.h>
#include <fmt/ranges.h>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/vector.h>
#include <armadillo>

#include "dense.h"
#include "gelex/data/bed_reader.h"
#include "gelex/data/grm.h"
#include "gelex/estimator/freq/estimator.h"
#include "gelex/model/freq/model.h"
#include "gelex/predictor.h"
#include "sparse.h"

namespace bind
{
namespace nb = nanobind;
using nb::literals::operator""_a;
namespace gx = gelex;

using arma::dmat;
using arma::dvec;

void gblup(nb::module_& m)
{
    nb::class_<gx::GBLUP>(m, "_GBLUP")
        .def(
            nb::init<std::string, dvec>(),
            "formula"_a,
            "phenotype"_a,
            nb::keep_alive<1, 3>())

        .def_prop_ro("n_individuals", &gx::GBLUP::n_individuals)

        .def_prop_ro("n_fixed_effects", &gx::GBLUP::n_fixed_effects)
        .def_prop_ro("n_random_effects", &gx::GBLUP::n_random_effects)
        .def_prop_ro("n_genetic_effects", &gx::GBLUP::n_genetic_effects)
        .def_prop_ro("n_gxe_effects", &gx::GBLUP::n_gxe_effects)

        .def_prop_ro("formula", &gx::GBLUP::formula)
        .def_prop_ro(
            "phenotype",
            &gx::GBLUP::phenotype,
            nb::rv_policy::reference_internal)
        .def(
            "add_fixed_effect",
            &gx::GBLUP::add_fixed_effect,
            nb::keep_alive<1, 4>())
        .def("add_random_effect", &gx::GBLUP::add_random_effect)
        .def("add_genetic_effect", &gx::GBLUP::add_genetic_effect)
        .def(
            "add_gxe_effect",
            &gx::GBLUP::add_gxe_effect,
            "name"_a,
            "design_mat_genetic"_a,
            "genetic_cov_mat"_a,
            "design_mat_env"_a.noconvert())
        .def("add_residual", &gx::GBLUP::add_residual)
        .def("clear", &gelex::GBLUP::clear, "reset the model")
        .def(
            "__repr__",
            [](const gx::GBLUP& self)
            {
                return fmt::format(
                    "<GBLUP object at {:p}: {:d} Individuals, {:d} Common"
                    "effects, {:d} Random effects, {:d} Genetic effects",
                    static_cast<const void*>(&self),
                    self.n_individuals(),
                    self.n_fixed_effects(),
                    self.n_random_effects(),
                    self.n_genetic_effects());
            })
        .def(
            "__str__",
            [](const gx::GBLUP& self)
            {
                return fmt::format(
                    "┌─ GBLUP Model ─────────────────────────────────\n"
                    "│ Individuals:    {:6d}\n"
                    "│ Common Effects:    {:6d}\n"
                    "│ Random Effects:    {:6d}\n"
                    "│ Genetic Effects:  {:6d}\n"
                    "└───────────────────────────────────────────────",
                    self.n_individuals(),
                    self.n_fixed_effects(),
                    self.n_random_effects(),
                    self.n_genetic_effects());
            });
}

void predictor(nb::module_& m)
{
    nb::class_<gx::Predictor>(m, "_Predictor")
        .def(
            nb::init<std::string_view, gelex::GBLUPParams>(),
            "train_bed"_a,
            "params"_a,
            nanobind::keep_alive<1, 3>())
        .def(
            "set_cross_grm",
            &gx::Predictor::set_cross_grm,
            nb::keep_alive<1, 3>())
        .def("_compute_random_effects", &gx::Predictor::compute_random_effects)
        .def("_compute_fixed_effects", &gx::Predictor::compute_fixed_effects)
        .def_prop_ro("test_individuals", &gelex::Predictor::test_individuals);
}

void estimator(nb::module_& m)
{
    nb::class_<gx::Estimator>(m, "Estimator")
        .def(
            nb::init<std::string_view, size_t, double>(),
            "optimizer"_a = "AI",
            "max_iter"_a = 20,
            "tol"_a = 1e-8,
            "Initialize the Estimator\n\n"
            "Parameters\n"
            "----------\n"
            "optimizer : str, optional\n"
            "    The optimization algorithm to use (default: 'NR')\n"
            "max_iter : int, optional\n"
            "    Maximum number of iterations (default: 20)\n"
            "tol : float, optional\n"
            "    Convergence tolerance (default: 1e-8)")
        .def(
            "fit",
            &gx::Estimator::fit,
            "model"_a,
            "em_init"_a = true,
            "verbose"_a = true,
            "fit the model\n\n"
            "Parameters\n"
            "----------\n"
            "model : GBLUP\n"
            "    The linear mixed model to fit\n"
            "em_init : bool, optional\n"
            "    Whether to use EM algorithm for initialization (default: "
            "True)\n\n"
            "verbose : bool, optional\n"
            "    Whether to print the optimization process (default: True)\n\n"
            "Returns\n"
            "-------\n"
            "None")
        .def(
            "set_optimizer",
            &gx::Estimator::set_optimizer,
            "optimizer"_a = "NR",
            "tol"_a = 1e-8,
            "reset optimizer\n\n"
            "Parameters\n"
            "----------\n"
            "optimizer : str, optional\n"
            "    The optimization algorithm to use (default: 'NR')\n"
            "max_iter : int, optional\n"
            "    Maximum number of iterations (default: 20)\n"
            "tol : float, optional\n"
            "    Convergence tolerance (default: 1e-8)\n\n"
            "Returns\n"
            "-------\n"
            "None");
}

void gblup_params(nb::module_& m)
{
    nb::class_<gx::GBLUPParams>(m, "GBLUPParams")
        .def(nb::init<dvec, dvec, dvec, std::vector<std::string>>());
}

void add_grm(nb::module_& m)
{
    nb::class_<gx::AddGrm>(m, "add_grm")
        .def(
            nb::init<std::string_view, size_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
            "exclude_individuals"_a = std::vector<std::string>{},
            "Additive Genomic Relationship Matrix calculation.\n\n"
            "Parameters\n"
            "----------\n"
            "bed_file: str\n"
            "    The plink bed file path\n"
            "chunk_size: int, optional\n"
            "    Number of snps to processed per step (default: 10000)\n\n"
            "Returns\n"
            "-------\n"
            "np.ndarray\n")
        .def("compute", &gx::AddGrm::compute, nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](gx::AddGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center", &gx::DomGrm::center, nb::rv_policy::reference_internal)
        .def_prop_ro("scale_factor", &gx::DomGrm::scale_factor);
}

void dom_grm(nb::module_& m)
{
    nb::class_<gx::DomGrm>(m, "dom_grm")
        .def(
            nb::init<std::string_view, size_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = gelex::DEFAULT_CHUNK_SIZE,
            "exclude_individuals"_a = std::vector<std::string>{},
            "Dominance Genomic Relationship Matrix calculation.\n\n"
            "Parameters\n"
            "----------\n"
            "bed_file: str\n"
            "    The plink bed file path\n"
            "chunk_size: int, optional\n"
            "    Number of snps to processed per step (default: 10000)\n\n"
            "Returns\n"
            "-------\n"
            "np.ndarray\n")
        .def("compute", &gx::DomGrm::compute, nb::rv_policy::move)
        .def_prop_ro(
            "individuals",
            [](gx::DomGrm& self) { return self.bed().individuals(); },
            nb::rv_policy::reference_internal)
        .def_prop_ro(
            "center", &gx::DomGrm::center, nb::rv_policy::reference_internal)
        .def_prop_ro("scale_factor", &gx::DomGrm::scale_factor);
}

void bed_reader(nb::module_& m)
{
    nb::class_<gx::BedReader>(m, "_BedReader")
        .def(
            nb::init<std::string_view, size_t, std::vector<std::string>>(),
            "bed_file"_a,
            "chunk_size"_a = gx::DEFAULT_CHUNK_SIZE,
            "dropped_ids"_a = std::vector<std::string>{},
            "Read a BED file in chunks.\n\n"
            "Parameters\n"
            "----------\n"
            "bed_file: str\n"
            "    The plink bed file path\n"
            "chunk_size: int, optional\n"
            "    Number of snps to processed per step (default: 10000)\n\n"
            "Returns\n"
            "-------\n"
            "None")
        .def("read_chunk", &gx::BedReader::read_chunk, nb::rv_policy::move)
        .def_prop_ro("num_snps", &gx::BedReader::num_snps)
        .def_prop_ro("num_individuals", &gx::BedReader::num_individuals)
        .def_prop_ro("snps", &gx::BedReader::snps)
        .def_prop_ro("individuals", &gx::BedReader::individuals);
}

void sp_dense_dot(nb::module_& m)
{
    m.def(
        "_sp_dense_dot",
        [](const arma::sp_dmat& a, const arma::dmat& b) -> arma::dmat
        {
            if (gelex::check_eye(a))
            {
                return b;
            }
            return a * b;
        },
        "a"_a,
        "b"_a,
        "Sparse Dense multiplication\n\n"
        "Parameters\n"
        "----------\n"
        "a : csc_matrix\n"
        "b : np.ndarray\n"
        "Returns\n"
        "-------\n"
        "np.ndarray\n"
        "    Resulting vector (individuals x 1)");
}
}  // namespace bind
