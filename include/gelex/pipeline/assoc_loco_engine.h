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

#ifndef GELEX_PIPELINE_ASSOC_LOCO_ENGINE_H_
#define GELEX_PIPELINE_ASSOC_LOCO_ENGINE_H_

#include "gelex/infra/logging/assoc_event.h"
#include "gelex/infra/logging/reml_event.h"
#include "gelex/pipeline/assoc_normal_engine.h"

namespace gelex
{

class PhenoPipe;
class GrmPipe;

class AssocLocoEngine
{
   public:
    using Config = AssocNormalEngine::Config;

    explicit AssocLocoEngine(Config config);

    auto run(
        PhenoPipe& pheno,
        GrmPipe& grm,
        const AssocObserver& observer = {},
        const RemlObserver& reml_observer = {}) -> void;

   private:
    Config config_;
};

}  // namespace gelex

#endif  // GELEX_PIPELINE_ASSOC_LOCO_ENGINE_H_
