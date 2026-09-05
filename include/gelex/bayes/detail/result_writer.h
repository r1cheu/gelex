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

#ifndef GELEX_BAYES_DETAIL_RESULT_WRITER_H_
#define GELEX_BAYES_DETAIL_RESULT_WRITER_H_

#include "gelex/bayes/genetic/gaussian.h"
#include "gelex/bayes/genetic/joint_spike_slab.h"
#include "gelex/bayes/genetic/scaled_mixture.h"
#include "gelex/bayes/genetic/spike_slab.h"
#include "gelex/bayes/mode_values.h"
#include "gelex/genetic_mode.h"
#include "gelex/io/detail/text_writer.h"

namespace gelex::detail
{

template <GeneticModeSet Modes, typename... Results>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const ModeValues<Modes, Results...>& result) -> void
{
    result.for_each([&]<GeneticMode Mode>(const auto& mode_result)
                    { write_family_summary_rows(writer, mode_result); });
}

template <typename ModeValuesType, typename JointResult>
auto write_genetic_summary_rows(
    TextWriter& writer,
    const JointModeValues<ModeValuesType, JointResult>& result) -> void
{
    write_genetic_summary_rows(writer, result.mode_values());
    write_family_summary_rows(writer, result.joint());
}

}  // namespace gelex::detail

#endif  // GELEX_BAYES_DETAIL_RESULT_WRITER_H_
