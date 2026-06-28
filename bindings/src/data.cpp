#include <nanobind/eigen/dense.h>
#include <nanobind/nanobind.h>

#include "gelex/data/genotype_method.h"
#include "gelex/data/locus_encoding.h"
#include "gelex/types/genetic_effect_type.h"

#include "gelex_py/register.h"

namespace nb = nanobind;

namespace gelex_py
{

void register_data(nb::module_& m)
{
    using namespace nb::literals;

    m.def(
        "encode_inplace",
        &gelex::encode_inplace<double>,
        "genotypes"_a,
        "effect"_a,
        "method"_a,
        "tol"_a = 1e-12,
        "marker_offset"_a = 0,
        "Encode a genotype matrix in place and return its LociEncoding.\n\n"
        "genotypes must be a writable float64 NumPy array in column-major\n"
        "(order='F') layout with shape (n_samples, n_markers); raw genotypes\n"
        "are 0/1/2 with NaN for missing. The array is modified in place (no\n"
        "copy); any other dtype/layout raises a TypeError.");
}

}  // namespace gelex_py
