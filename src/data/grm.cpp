#include "gelex/data/grm.h"

#include <memory>

#include <Eigen/Core>
#ifdef USE_MKL
#include <mkl.h>
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif

#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{
using Eigen::Index;

GRM::GRM(std::filesystem::path bed_path)
    : bed_(
          bed_path,
          std::make_shared<SampleManager>(
              bed_path.replace_extension(".fam"),
              false))
{
}

auto GRM::update_grm(
    Eigen::Ref<Eigen::MatrixXd> grm,
    const Eigen::Ref<const Eigen::MatrixXd>& genotype) -> void
{
    const auto n = static_cast<int>(genotype.rows());
    const auto m = static_cast<int>(genotype.cols());

    // dsyrk: C := alpha * A * A^T + beta * C
    // CblasLower: only updates lower triangle
    cblas_dsyrk(
        CblasColMajor,
        CblasLower,
        CblasNoTrans,
        n,
        m,
        1.0,
        genotype.data(),
        n,
        1.0,
        grm.data(),
        n);
}
}  // namespace gelex
