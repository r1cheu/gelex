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

#include "predict_pipe.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "gelex/data/bed_pipe.h"
#include "gelex/data/sample_manager.h"

namespace gelex
{

PredictDataPipe::PredictDataPipe(const Config& config)
{
    auto fam_path = config.bed_path;
    fam_path.replace_extension(".fam");
    sample_manager_ = std::make_shared<SampleManager>(fam_path);

    if (!config.qcovar_path.empty())
    {
        load_qcovariates(config);
    }

    if (!config.dcovar_path.empty())
    {
        load_dcovariates(config);
    }

    intersect();
    format_dcovariates();

    load_genotype(config);
}

auto PredictDataPipe::load_qcovariates(const Config& config) -> void
{
    qcovar_frame_ = DataFrame<double>::read(config.qcovar_path);
    qcovariate_names_ = qcovar_frame_->columns();
}

auto PredictDataPipe::load_dcovariates(const Config& config) -> void
{
    dcovar_frame_ = DataFrame<std::string>::read(config.dcovar_path);
    dcovariate_names_ = dcovar_frame_->columns();
}

auto PredictDataPipe::load_genotype(const Config& config) -> void
{
    BedPipe bed_pipe(config.bed_path, sample_manager_);
    genotypes_ = bed_pipe.load();
}

auto PredictDataPipe::intersect() -> void
{
    if (qcovar_frame_)
    {
        sample_manager_->intersect(qcovar_frame_->index_column().data());
    }

    if (dcovar_frame_)
    {
        sample_manager_->intersect(dcovar_frame_->index_column().data());
    }

    sample_manager_->finalize();
}

auto PredictDataPipe::format_dcovariates() -> void
{
    const auto num_samples
        = static_cast<Eigen::Index>(sample_manager_->num_common_samples());
    const auto& common_ids = sample_manager_->common_ids();
    sample_ids_ = common_ids;

    if (qcovar_frame_)
    {
        auto aligned = *qcovar_frame_;
        aligned.intersect_index_inplace(common_ids);
        auto qcov = aligned.eigen();

        qcovariates_ = Eigen::MatrixXd::Ones(num_samples, qcov.cols() + 1);
        if (qcov.size() != 0)
        {
            qcovariates_.middleCols(1, qcov.cols()) = qcov;
        }
    }
    else
    {
        qcovariates_ = Eigen::MatrixXd::Ones(num_samples, 1);
    }

    if (dcovar_frame_)
    {
        auto aligned = *dcovar_frame_;
        aligned.intersect_index_inplace(common_ids);

        dcovariates_.clear();
        const auto names = aligned.columns();
        for (size_t col = 0; col < names.size(); ++col)
        {
            dcovariates_[names[col]] = aligned.column(col).data();
        }
    }
}

}  // namespace gelex
