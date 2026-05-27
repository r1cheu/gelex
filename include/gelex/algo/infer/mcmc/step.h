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

#ifndef GELEX_ALGO_INFER_MCMC_STEP_H_
#define GELEX_ALGO_INFER_MCMC_STEP_H_

namespace gelex::mcmc
{

class Step
{
   public:
    Step() = default;
    Step(const Step&) = delete;
    Step(Step&&) noexcept = delete;
    auto operator=(const Step&) -> Step& = delete;
    auto operator=(Step&&) noexcept -> Step& = delete;
    virtual ~Step() = default;

    virtual auto step() -> void = 0;
};

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_STEP_H_
