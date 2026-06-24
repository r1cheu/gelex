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

#ifndef APPS_CLI_PREDICT_REPORTER_H_
#define APPS_CLI_PREDICT_REPORTER_H_

#include <cstddef>
#include <string_view>

#include "gelex/data/genotype_method.h"

namespace cli
{

class PredictReporter
{
   public:
    auto show_snp_selection(
        size_t num_matched,
        size_t num_missing,
        size_t num_mismatched,
        size_t num_total,
        std::string_view bfile_path,
        std::string_view snp_effect_path) const -> void;
    auto show_data_loaded(
        size_t num_samples,
        size_t num_snps,
        size_t num_covar_terms,
        gelex::GenotypeMethod geno_method) const -> void;
    auto show_results_written(std::string_view output_path, size_t num_samples)
        const -> void;
};

}  // namespace cli

#endif  // APPS_CLI_PREDICT_REPORTER_H_
