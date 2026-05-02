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

#ifndef GELEX_ALGO_INFER_VI_CHAIN_H_
#define GELEX_ALGO_INFER_VI_CHAIN_H_

#include <tuple>
#include <utility>

namespace gelex::vi
{

template <typename... Samplers>
class CaviChain
{
   public:
    static_assert(
        sizeof...(Samplers) > 0,
        "CaviChain must contain at least one sampler");

    explicit CaviChain(Samplers... samplers) : samplers_(std::move(samplers)...)
    {
    }

    CaviChain(const CaviChain&) = delete;
    auto operator=(const CaviChain&) -> CaviChain& = delete;
    CaviChain(CaviChain&&) noexcept = default;
    auto operator=(CaviChain&&) -> CaviChain& = delete;
    ~CaviChain() = default;

    auto step() -> void
    {
        std::apply(
            [](auto&... sampler) { (sampler.sample(), ...); }, samplers_);
    }

   private:
    std::tuple<Samplers...> samplers_;
};

template <typename... Samplers>
CaviChain(Samplers...) -> CaviChain<Samplers...>;

}  // namespace gelex::vi

#endif  // GELEX_ALGO_INFER_VI_CHAIN_H_
