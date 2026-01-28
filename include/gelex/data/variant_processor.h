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

#ifndef GELEX_DATA_VARIANT_PROCESSOR_H_
#define GELEX_DATA_VARIANT_PROCESSOR_H_

#include <concepts>

#include <Eigen/Dense>

namespace gelex
{

struct VariantStats
{
    double mean{0.0};
    double stddev{0.0};
    bool is_monomorphic{false};
};

template <typename T>
concept VariantProcessor
    = requires(T processor, Eigen::Ref<Eigen::VectorXd> variant) {
          { T::process_variant(variant) } -> std::same_as<VariantStats>;
      };

namespace detail
{
struct GenotypeCounts
{
    double p_AA{0.0};
    double p_Aa{0.0};
    double p_aa{0.0};
};

VariantStats compute_and_standardize(Eigen::Ref<Eigen::VectorXd> variant);
GenotypeCounts count_frequencies(
    const Eigen::Ref<const Eigen::VectorXd>& variant);

}  // namespace detail
}  // namespace gelex

namespace gelex
{
struct StandardizingProcessor
{
    static constexpr bool dom = false;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct RawProcessor
{
    static constexpr bool dom = false;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct HardWenbergProcessor
{
    static constexpr bool dom = false;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct NOIAProcessor
{
    static constexpr bool dom = false;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantStandardizingProcessor
{
    static constexpr bool dom = true;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantRawProcessor
{
    static constexpr bool dom = true;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantOrthogonalHWEProcessor
{
    static constexpr bool dom = true;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

struct DominantNOIAProcessor
{
    static constexpr bool dom = true;
    static VariantStats process_variant(Eigen::Ref<Eigen::VectorXd> variant);
};

}  // namespace gelex

#endif  // GELEX_DATA_VARIANT_PROCESSOR_H_
