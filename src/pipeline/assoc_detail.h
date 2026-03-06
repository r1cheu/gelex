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

#ifndef GELEX_PIPELINE_ASSOC_DETAIL_H_
#define GELEX_PIPELINE_ASSOC_DETAIL_H_

#include <cstddef>
#include <functional>

#include <Eigen/Core>

#include "gelex/data/genotype/bed_pipe.h"
#include "gelex/infra/logging/assoc_event.h"
#include "gelex/types/assoc_input.h"
#include "gelex/types/chr_group.h"

namespace gelex
{
class FreqModel;
class FreqState;
enum class GenotypeProcessMethod : uint8_t;
enum class ModelType : uint8_t;

namespace detail
{

auto dispatch_assoc_chunk_by_method(
    GenotypeProcessMethod method,
    ModelType model_type,
    Eigen::Ref<Eigen::MatrixXd> genotype,
    Eigen::VectorXd* freqs) -> void;

auto update_assoc_input(
    AssocInput& input,
    const FreqModel& model,
    const FreqState& state,
    Eigen::MatrixXd&& v_inv) -> void;

class ChrScanner
{
   public:
    struct Config
    {
        int chunk_size;
        size_t total_snps;
    };

    struct SnpResult
    {
        double freq;
        double beta;
        double se;
        double p_value;
    };

    using GenoProcessor
        = std::function<void(Eigen::Ref<Eigen::MatrixXd>, Eigen::VectorXd*)>;
    using ResultWriter
        = std::function<void(size_t snp_index, const SnpResult&)>;

    ChrScanner(Config config, BedPipe& bed, AssocObserver observer);

    ChrScanner(const ChrScanner&) = delete;
    auto operator=(const ChrScanner&) -> ChrScanner& = delete;

    auto scan(
        const ChrGroup& group,
        const GenoProcessor& processor,
        const ResultWriter& writer) -> void;

    auto assoc_input() -> AssocInput& { return input_; }

   private:
    Config config_;
    BedPipe& bed_;
    AssocObserver observer_;

    AssocInput input_;
    AssocOutput output_;
    Eigen::VectorXd freqs_;
    size_t progress_counter_ = 0;
};

}  // namespace detail
}  // namespace gelex

#endif  // GELEX_PIPELINE_ASSOC_DETAIL_H_
