/*
 * Copyright 2026 RuLei Chen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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

namespace
{
auto create_sample_manager(const std::filesystem::path& bed_path)
    -> std::shared_ptr<SampleManager>
{
    auto fam_path = std::filesystem::path(bed_path).replace_extension(".fam");
    auto manager = std::make_shared<SampleManager>(fam_path, false);
    manager->finalize();
    return manager;
}
}  // namespace

GRM::GRM(std::filesystem::path bed_path)
    : sample_manager_(create_sample_manager(bed_path)),
      bed_(bed_path, sample_manager_)
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
