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

#ifndef GELEX_ALGO_INFER_MCMC_CHAIN_H_
#define GELEX_ALGO_INFER_MCMC_CHAIN_H_

#include <tuple>
#include <utility>

namespace gelex::mcmc
{

template <typename... Samplers>
class McmcChain
{
   public:
    static_assert(
        sizeof...(Samplers) > 0,
        "McmcChain must contain at least one sampler");

    explicit McmcChain(Samplers... samplers) : samplers_(std::move(samplers)...)
    {
    }

    McmcChain(const McmcChain&) = delete;
    auto operator=(const McmcChain&) -> McmcChain& = delete;
    McmcChain(McmcChain&&) noexcept = default;
    auto operator=(McmcChain&&) -> McmcChain& = delete;
    ~McmcChain() = default;

    auto step() -> void
    {
        std::apply(
            [](auto&... sampler) { (sampler.sample(), ...); }, samplers_);
    }

   private:
    std::tuple<Samplers...> samplers_;
};

template <typename... Samplers>
McmcChain(Samplers...) -> McmcChain<Samplers...>;

}  // namespace gelex::mcmc

#endif  // GELEX_ALGO_INFER_MCMC_CHAIN_H_
