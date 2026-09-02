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

#include <gelex/bayes/stats/diagnostics.h>
#include <gelex/genetic_mode.h>

auto main() -> int
{
    const auto modes = gelex::GeneticMode::A | gelex::GeneticMode::D;
    return modes.size() == 2 && gelex::fft_next_fast_len(7) == 8 ? 0 : 1;
}
