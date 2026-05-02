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

#ifndef GELEX_MODEL_BAYES_CHECKPOINT_READER_H_
#define GELEX_MODEL_BAYES_CHECKPOINT_READER_H_

#include <filesystem>

#include "gelex/model/bayes/checkpoint.h"

namespace gelex
{

[[nodiscard]] auto read_checkpoint(const std::filesystem::path& path)
    -> Checkpoint;

}  // namespace gelex

#endif  // GELEX_MODEL_BAYES_CHECKPOINT_READER_H_
