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

#ifndef GELEX_ALGO_INFER_CHAIN_H_
#define GELEX_ALGO_INFER_CHAIN_H_

#include <tuple>
#include <utility>

namespace gelex::infer
{

template <typename... Steps>
class Chain
{
   public:
    static_assert(sizeof...(Steps) > 0, "Chain must contain at least one step");

    explicit Chain(Steps... steps) : steps_(std::move(steps)...) {}

    Chain(const Chain&) = delete;
    auto operator=(const Chain&) -> Chain& = delete;
    Chain(Chain&&) noexcept = default;
    auto operator=(Chain&&) -> Chain& = delete;
    ~Chain() = default;

    auto step() -> void
    {
        std::apply([](auto&... s) { (s.step(), ...); }, steps_);
    }

   private:
    std::tuple<Steps...> steps_;
};

template <typename... Steps>
Chain(Steps...) -> Chain<Steps...>;

}  // namespace gelex::infer

#endif  // GELEX_ALGO_INFER_CHAIN_H_
